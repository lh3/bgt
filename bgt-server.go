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
		_, present := r.Form["G"];
		if present {
			flag |= 2; // BGT_F_NO_GT
		}
		_, present = r.Form["C"];
		if present {
			flag |= 1; // BGT_F_SET_AC
	  	}
		_, present = r.Form["S"];
		if present {
			flag |= 4; // BGT_F_CNT_HAP
	  	}
		_, present = r.Form["H"];
		if present {
			flag |= 8; // BGT_F_CNT_HAP
	  	}
		_, present = r.Form["s"];
		if present {
			flag |= 1; // BGT_F_SET_AC
		}
		C.bgtm_set_flag(bm, C.int(flag));
		if (flag & 12) != 0 {
			vcf_out = false;
		}
	}
	{ // set site filter
		a, present := r.Form["f"];
		if present {
			cstr := C.CString(a[0]);
			C.bgtm_set_flt_site(bm, cstr);
			C.free(unsafe.Pointer(cstr));
		}
	}
	{ // set region
		a, present := r.Form["r"];
		if present {
			cstr := C.CString(a[0]);
			C.bgtm_set_region(bm, cstr);
			C.free(unsafe.Pointer(cstr));
		}
	}
	{ // set start
		a, present := r.Form["i"];
		if present {
			i, _ := strconv.Atoi(a[0]);
			C.bgtm_set_start(bm, C.int64_t(i));
		}
	}
	{ // set start
		a, present := r.Form["n"];
		if present {
			max_read, _ = strconv.Atoi(a[0]);
		}
	}
	{ // set alleles
		a, present := r.Form["a"];
		if present {
			cstr := C.CString(a[0]);
			C.bgtm_set_alleles(bm, cstr, nil, nil);
			C.free(unsafe.Pointer(cstr));
		}
	}
	{ // set sample groups
		a, present := r.Form["s"];
		if present {
			for _, s := range a {
				cstr := C.CString(s);
				C.bgtm_add_group(bm, cstr);
				C.free(unsafe.Pointer(cstr));
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
		if n_read >= max_read {
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
		}
		n_read += 1;
	}
	C.bcf_destroy1(b);

	// print hapcnt
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
}

/*****************
 * Main function *
 *****************/

func main() {
	port := ":8000";

	// parse command line options
	for {
		opt, arg := getopt(os.Args, "p:");
		if opt == 'p' {
			port = fmt.Sprintf(":%s", arg);
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

	http.HandleFunc("/query", bgs_query);
	fmt.Fprintln(os.Stderr, "MESSAGE: server started...");
	http.ListenAndServe(port, nil);
}
