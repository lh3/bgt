package main

import (
	"os"
	"fmt"
	"net/http"
	"unsafe"
	"strconv"
	"strings"
)

/*
#cgo LDFLAGS: -L. -lbgt -lpthread -lz -lm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bgt.h"

bgtm_t *bgtm_reader_init_n(int n_files)
{
	bgtm_t *bm;
	bm = (bgtm_t*)calloc(1, sizeof(bgtm_t));
	bm->n_bgt = n_files;
	bm->bgt = (bgt_t**)calloc(bm->n_bgt, sizeof(void*));
	bm->r = (bgt_rec_t*)calloc(bm->n_bgt, sizeof(bgt_rec_t));
	return bm;
}

void bgtm_set_file(bgtm_t *bm, int i, const bgt_file_t *f)
{
	bm->bgt[i] = bgt_reader_init(f);
}

char *bgtm_format_bcf1(const bgtm_t *bm, const bcf1_t *v)
{
	kstring_t s = {0,0,0};
	vcf_format1(bm->h_out, v, &s);
	return s.s;
}

char *bgtm_hapcnt2str(const bgtm_t *bm)
{
	bgt_hapcnt_t *hc;
	int n_hap;
	char *s;
	hc = bgtm_hapcnt(bm, &n_hap);
	return bgtm_hapcnt_print_destroy(bm, n_hap, hc);
}
*/
import "C"

/****************
 * BSD getopt() *
 ****************/

var optind int = 1
var getopt_place int = -1

func getopt(args []string, ostr string) (int, string) {
	if getopt_place == -1 { // update scanning pointer
		if optind >= len(args) || args[optind][0] != '-' {
			getopt_place = -1
			return -1, ""
		}
		if optind < len(args) {
			getopt_place = 0
		}
		if getopt_place + 1 < len(args[optind]) {
			getopt_place += 1
			if args[optind][getopt_place] == '-' { // found "--"
				optind += 1
				getopt_place = -1
				return -1, ""
			}
		}
	}
	optopt := args[optind][getopt_place];
	getopt_place += 1
	oli, arg := strings.IndexByte(ostr, optopt), "";
	if optopt == ':' || oli < 0 {
		if optopt == '-' {
			return -1, ""
		}
		if getopt_place < 0 {
			optind += 1
		}
		return '?', ""
	}
	if oli + 1 >= len(ostr) || ostr[oli+1] != ':' {
		if getopt_place < 0 || getopt_place >= len(args[optind]) {
			optind += 1
			getopt_place = -1
		}
	} else {
		if getopt_place >= 0 && getopt_place < len(args[optind]) {
			arg = args[optind][getopt_place:]
		} else if optind += 1; len(args) <= optind { // no arguments
			getopt_place = -1
			if len(ostr) > 0 && ostr[0] == ':' {
				return ':', ""
			}
			return '?', ""
		} else {
			arg = args[optind]
		}
		getopt_place = -1
		optind += 1
	}
	return int(optopt), arg
}

/****************
 * BGT wrappers *
 ****************/

var bgt_files [](*C.bgt_file_t);
var bgt_vardb *C.fmf_t = nil;
var bgt_file_names []string;
var bgt_port string = ":8000";
var bgt_max_gt uint64 = uint64(1000000);

func bgtm_open(fns []string) ([](*C.bgt_file_t)) {
	files := make([](*C.bgt_file_t), len(fns));
	for i := 0; i < len(fns); i += 1 {
		cstr := C.CString(fns[i]);
		defer C.free(unsafe.Pointer(cstr));
		files[i] = C.bgt_open(cstr);
	}
	return files;
}

func bgtm_close(files [](*C.bgt_file_t)) {
	for i := 0; i < len(files); i += 1 {
		C.bgt_close(files[i]);
	}
}

func bgtm_reader_init(files [](*C.bgt_file_t)) (*C.bgtm_t) {
	bm := C.bgtm_reader_init_n(C.int(len(files)));
	for i := 0; i < len(files); i += 1 {
		C.bgtm_set_file(bm, C.int(i), files[i]);
	}
	return bm;
}

/************
 * Handlers *
 ************/

func bgs_replace_op(t string) string {
	s := strings.Replace(t, ".AND.", "&&", -1);
	s = strings.Replace(s, ".and.", "&&", -1);
	s = strings.Replace(s, ".OR.", "||", -1);
	s = strings.Replace(s, ".or.", "||", -1);
	return s;
}

func bgs_query(w http.ResponseWriter, r *http.Request) {
	r.URL.RawQuery = strings.Replace(r.URL.RawQuery, "&&", ".AND.", -1);
	r.ParseForm();
	flag := 0;
	max_read := 2147483647;
	vcf_out := true;
	bm := bgtm_reader_init(bgt_files);

	{ // set flag
		if len(r.Form["G"]) > 0 {
			flag |= 2; // BGT_F_NO_GT
		}
		if len(r.Form["C"]) > 0 || len(r.Form["s"]) > 0 {
			flag |= 1; // BGT_F_SET_AC
	  	}
		if len(r.Form["S"]) > 0 {
			flag |= 4; // BGT_F_CNT_HAP
	  	}
		if len(r.Form["H"]) > 0 {
			flag |= 8; // BGT_F_CNT_HAP
	  	}
		C.bgtm_set_flag(bm, C.int(flag));
		if (flag & 12) != 0 {
			vcf_out = false;
		}
	}
	if len(r.Form["f"]) > 0 { // set site filter
		cstr := C.CString(bgs_replace_op(r.Form["f"][0]));
		ret := int(C.bgtm_set_flt_site(bm, cstr));
		C.free(unsafe.Pointer(cstr));
		if ret != 0 {
			C.bgtm_reader_destroy(bm);
			http.Error(w, "400 Bad Request: failed to parse parameter 'f'", 400);
			return;
		}
	}
	if len(r.Form["r"]) > 0 { // set region
		cstr := C.CString(r.Form["r"][0]);
		ret := int(C.bgtm_set_region(bm, cstr));
		C.free(unsafe.Pointer(cstr));
		if ret < 0 {
			C.bgtm_reader_destroy(bm);
			http.Error(w, "400 Bad Request: failed to set region with parameter 'r'", 400);
			return;
		}
	}
	if len(r.Form["i"]) > 0 { // set start
		i, _ := strconv.Atoi(r.Form["i"][0]);
		C.bgtm_set_start(bm, C.int64_t(i));
	}
	if len(r.Form["n"]) > 0 { // set max number of records to read
		max_read, _ = strconv.Atoi(r.Form["n"][0]);
	}
	if len(r.Form["t"]) > 0 { // tabular output
		cstr := C.CString(r.Form["t"][0]);
		C.bgtm_set_table(bm, cstr);
		C.free(unsafe.Pointer(cstr));
		vcf_out = false;
	}
	if len(r.Form["a"]) > 0 { // set alleles
		cstr := C.CString(r.Form["a"][0]);
		n_al := int(C.bgtm_set_alleles(bm, cstr, bgt_vardb, nil));
		C.free(unsafe.Pointer(cstr));
		if n_al <= 0 {
			C.bgtm_reader_destroy(bm);
			if n_al < 0 {
				http.Error(w, "400 Bad Request: failed to retrieve alleles with parameter 'a'", 400);
			} else {
				http.Error(w, "204 No Content: no alleles matching parameter 'a'", 204);
			}
			return;
		}
	}
	if len(r.Form["s"]) > 0 { // set sample groups
		for _, s := range r.Form["s"] {
			cstr := C.CString(bgs_replace_op(s));
			ret := C.bgtm_add_group(bm, cstr);
			C.free(unsafe.Pointer(cstr));
			if ret < 0 {
				http.Error(w, "400 Bad Request: failed to set sample group with parameter 's'", 400);
				return;
			}
		}
	}
	C.bgtm_prepare(bm);

	// print header if necessary
	if vcf_out {
		gstr := C.GoString(bm.h_out.text);
		fmt.Fprintln(w, gstr);
	}

	// read through
	b := C.bcf_init1();
	n_read := 0;
	for {
		if n_read > max_read || uint64(bm.n_gt_read) > bgt_max_gt {
			break;
		}
		ret := int(C.bgtm_read(bm, b));
		if ret < 0 {
			break;
		}
		if vcf_out {
			s := C.bgtm_format_bcf1(bm, b);
			gstr := C.GoString(s);
			fmt.Fprintln(w, gstr);
			C.free(unsafe.Pointer(s));
		} else if int(bm.n_fields) > 0 {
			gstr := C.GoString(bm.tbl_line.s);
			fmt.Fprintln(w, gstr);
		}
		n_read += 1;
	}
	C.bcf_destroy1(b);

	// print hapcnt and/or sample list
	if !vcf_out && int(bm.n_aal) > 0 {
		if (flag & 8) != 0 {
			s := C.bgtm_hapcnt2str(bm);
			gstr := C.GoString(s);
			fmt.Fprint(w, gstr);
			C.free(unsafe.Pointer(s));
		}
		if (flag & 4) != 0 {
			s := C.bgtm_alcnt_print(bm);
			gstr := C.GoString(s);
			fmt.Fprint(w, gstr);
			C.free(unsafe.Pointer(s));
		}
	}

	C.bgtm_reader_destroy(bm);

	if n_read > max_read {
		fmt.Fprintln(w, "*");
	}
}

func bgs_help(w http.ResponseWriter, r *http.Request) {
	fmt.Fprintln(w, "Manual");
	fmt.Fprintln(w, "======\n");
	fmt.Fprintln(w, "The BGT server has a similar interface to 'bgt view'. Query examples:\n");
	fmt.Fprintf(w,  "  curl --data 's=(population==\"FIN\")&s=(population==\"CEU\")&G&f=(AC1>0.and.AC2>0)' 0.0.0.0%s/query\n", bgt_port);
	fmt.Fprintf(w,  "  curl --data 'a=,11:151344:1:G,11:110992:AACTT:A,11:160513::G&S&s=(population==\"FIN\")' 0.0.0.0%s/query\n\n", bgt_port);
	fmt.Fprintln(w, "Sample selection parameter:\n");
	fmt.Fprintln(w, "  s EXPR  List of samples in a comma-leading comma-separate list (e.g. ,sample1,sample2) or an");
	fmt.Fprintln(w, "          expression (e.g. s=population==\"FIN\"). There can be multiple 's' parameters. Each of");
	fmt.Fprintln(w, "          them defines a sample group.\n");
	fmt.Fprintln(w, "Site selection parameters:\n");
	fmt.Fprintln(w, "  r STR   Region in a format like '11:200,000-300,000'\n");
	fmt.Fprintln(w, "  i INT   Start from the i-th record; INT>0\n");
	fmt.Fprintln(w, "  n INT   Read at most INT records\n");
	fmt.Fprintln(w, "  a EXPR  List of alleles in a format similar to parameter 's'. An allele is specified by");
	fmt.Fprintln(w, "          chr:1basedPos:refLen:alleleSeq. Conditions may not work unless the server is launched with");
	fmt.Fprintln(w, "          a variant annotation database.\n");
	fmt.Fprintln(w, "  f EXPR  Filters on per sample group allele counts. EXPR could include AC (primary allele count),");
	fmt.Fprintln(w, "          AN (total called alleles), AC# (primary allele count of the #-th sample group) and AN#.\n");
	fmt.Fprintln(w, "VCF output parameters:\n");
	fmt.Fprintln(w, "  G       Don't output sample genotypes\n");
	fmt.Fprintln(w, "  C       Output AC and AN VCF INFO fields. This parameter is automatically set if 's' is applied.\n");
	fmt.Fprintln(w, "Non-VCF output parameters:\n");
	fmt.Fprintln(w, "  S       Output samples having requested alleles (requiring parameter 'a')\n");
	fmt.Fprintln(w, "  H       Output counts of haplotypes across requested alleles (requiring parameter 'a')\n");
	fmt.Fprintln(w, "  t STR   Comma-separated list of fields in tabular output. Accepted variables:");
	fmt.Fprintln(w, "          CHROM, POS, END, REF, ALT, AC, AN, AC#, AN# (# for a group number)\n");
}

/*****************
 * Main function *
 *****************/

func main() {
	// parse command line options
	for {
		opt, arg := getopt(os.Args, "d:p:m:");
		if opt == 'p' {
			bgt_port = fmt.Sprintf(":%s", arg);
		} else if opt == 'm' {
			bgt_max_gt, _ = strconv.ParseUint(arg, 10, 64);
		} else if opt == 'd' {
			cstr := C.CString(arg);
			bgt_vardb = C.fmf_read(cstr);
			C.free(unsafe.Pointer(cstr));
		} else if opt < 0 {
			break;
		}
	}
	if optind == len(os.Args) {
		fmt.Fprintln(os.Stderr, "Usage: bgt-server [options] <bgt.pre1> [...]");
		os.Exit(1);
	}

	bgt_files = bgtm_open(os.Args[optind:]);
	defer bgtm_close(bgt_files);

	http.HandleFunc("/", bgs_help);
	http.HandleFunc("/query", bgs_query);
	fmt.Fprintln(os.Stderr, "MESSAGE: server started...");
	http.ListenAndServe(bgt_port, nil);
}
