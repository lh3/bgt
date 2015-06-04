package main

import (
	"os"
	"fmt"
	"net/http"
	"unsafe"
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

bgt_t *bgtm_get_bgt(bgtm_t *bm, int i)
{
	return bm->bgt[i];
}

void bgtm_al_set_rid(bgtm_t *bm, bgt_allele_t *a)
{
	if (bm->h_out == 0) bgtm_prepare(bm);
	a->rid = bcf_name2id(bm->h_out, a->chr.s);
}

const fmf1_t *fmf_get_row(const fmf_t *f, int r)
{
	return &f->rows[r];
}
*/
import "C"

var bgt_files [](*C.bgt_file_t);

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

func bgtm_open(fns []string) ([](*C.bgt_file_t)) {
	files := make([](*C.bgt_file_t), len(fns));
	for i := 0; i < len(fns); i += 1 {
		cstr := C.CString(fns[i]);
		defer C.free(unsafe.Pointer(cstr)); // on Mac, this triggers a weird error
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

func bgs_getvcf(w http.ResponseWriter, r *http.Request) {
	bm := bgtm_reader_init(bgt_files);
	r.ParseForm();
	C.bgtm_reader_destroy(bm);
}

func bgs_getbgt(w http.ResponseWriter, r *http.Request) {
	for i := 0; i < len(bgt_files); i += 1 {
		gstr := C.GoString(bgt_files[i].prefix);
		fmt.Fprintf(w, "%d\t%s\n", i + 1, gstr);
	}
}

func bgs_getspl(w http.ResponseWriter, r *http.Request) {
	r.URL.RawQuery = strings.Replace(r.URL.RawQuery, "&&", ".AND.", -1);
	r.ParseForm();
	var expr *C.char = nil;
	_, present_se := r.Form["se"];
	if present_se {
		s := strings.Replace(r.Form["se"][0], ".AND.", "&&", -1);
		s = strings.Replace(s, ".and.", "&&", -1);
		s = strings.Replace(s, ".OR.", "||", -1);
		s = strings.Replace(s, ".or.", "&&", -1);
		expr = C.CString(s);
		defer C.free(unsafe.Pointer(expr));
	}

	var cstr *C.char = nil;
	_, present_a := r.Form["al"];
	if present_a {
		al := make([]C.bgt_allele_t, len(r.Form["al"]));
		min_pos := 1<<30;
		max_pos := 0;
		for i := 0; i < len(al); i += 1 { // parse al
			cstr = C.CString(r.Form["al"][i]);
			C.bgt_al_parse(cstr, &al[i]);
			C.free(unsafe.Pointer(cstr));
			if i > 0 && C.strcmp(al[i].chr.s, al[0].chr.s) != 0 { // TODO: for now, only working when all al on the same chr
				i -= 1;
			}
			if min_pos > int(al[i].pos) {
				min_pos = int(al[i].pos);
			}
			if max_pos < int(al[i].pos) {
				max_pos = int(al[i].pos);
			}
		}
		bm := bgtm_reader_init(bgt_files);

		reg := fmt.Sprintf("%s:%d-%d", C.GoString(al[0].chr.s), min_pos+1, max_pos+1);
		cstr = C.CString(reg);
		C.bgtm_set_region(bm, cstr);
		C.free(unsafe.Pointer(cstr));

		if present_se {
			C.bgtm_add_group_core(bm, 0, nil, expr);
		}

		C.bgtm_prepare(bm);
		for i := 0; i < len(al); i += 1 { // parse al
			C.bgtm_al_set_rid(bm, &al[i]);
		}
		b := C.bcf_init1();
		match_cnt := make([]int, int(bm.n_out)>>1);
		for {
			ret := int(C.bgtm_read(bm, b));
			if ret < 0 {
				break;
			}
			ret = 0;
			for i := 0; i < len(al); i += 1 {
				ret = int(C.bgt_al_test(b, &al[i]));
				if ret != 0 {
					break;
				}
			}
			if ret == 0 {
				continue;
			}
			is_ref := (ret == 2);
			a0 := (*[1<<30]C.uint8_t)(unsafe.Pointer(bm.a[0]));
			a1 := (*[1<<30]C.uint8_t)(unsafe.Pointer(bm.a[1]));
			for i := 0; i < int(bm.n_out)>>1; i += 1 {
				g1 := a0[i<<1|0] | a1[i<<1|0]<<1;
				g2 := a0[i<<1|1] | a1[i<<1|1]<<1;
				if (is_ref) {
					if g1 == 0 || g2 == 0 {
						match_cnt[i] += 1;
					}
				} else {
					if g1 == 1 || g2 == 1 {
						match_cnt[i] += 1;
					}
				}
			}
		}
		C.bcf_destroy1(b);

		sidx := (*[1<<30]C.uint64_t)(unsafe.Pointer(bm.sample_idx));
		for i := 0; i < int(bm.n_out)>>1; i += 1 {
			if match_cnt[i] == len(al) {
				bgt := C.bgtm_get_bgt(bm, C.int(int64(sidx[i<<1])>>32));
				gstr := C.GoString(C.fmf_get_row(bgt.f.f, C.int(sidx[i<<1])).name);
				fmt.Fprintf(w, "%s\t%d\n", gstr, int(sidx[i<<1]>>32) + 1);
			}
		}

		C.bgtm_reader_destroy(bm);
		for i := 0; i < len(al); i += 1 { // free al
			C.free(unsafe.Pointer(al[i].chr.s));
		}
	} else {
		var ke *C.kexpr_t = nil;
		if present_se {
			var err C.int;
			ke = C.ke_parse(expr, &err);
		}
		for i := 0; i < len(bgt_files); i += 1 {
			f := bgt_files[i].f;
			for j := 0; j < int(f.n_rows); j += 1 {
				if ke == nil || C.fmf_test(f, C.int(j), ke) != C.int(0) {
					gstr := C.GoString(C.fmf_get_row(f, C.int(j)).name);
					fmt.Fprintf(w, "%s\t%d\n", gstr, i + 1);
				}
			}
		}
		if ke != nil {
			C.ke_destroy(ke);
		}
	}
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
		fmt.Println("Usage: bgt-server [options] <bgt.pre1> [...]");
		os.Exit(1);
	}

	bgt_files = bgtm_open(os.Args[optind:]);
	defer bgtm_close(bgt_files);

	http.HandleFunc("/getvcf", bgs_getvcf);
	http.HandleFunc("/getbgt", bgs_getbgt);
	http.HandleFunc("/getsamples", bgs_getspl);
	http.ListenAndServe(port, nil);
}
