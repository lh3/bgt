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
	fmt.Fprintf(w, "s=%s\n", r.Form["s"]);
	fmt.Fprintf(w, "s=%d\n", len(r.Form["s"]));
	C.bgtm_reader_destroy(bm);
}

func bgs_getbgt(w http.ResponseWriter, r *http.Request) {
	for i := 0; i < len(bgt_files); i += 1 {
		gstr := C.GoString(bgt_files[i].prefix);
		fmt.Fprintf(w, "%d\t%s\n", i + 1, gstr);
	}
}

func bgs_getspl(w http.ResponseWriter, r *http.Request) {
	r.ParseForm();
	_, present := r.Form["q"];
	if !present {
		for i := 0; i < len(bgt_files); i += 1 {
			f := bgt_files[i].f;
			for j := 0; j < int(f.n_rows); j += 1 {
				gstr := C.GoString(C.fmf_get_row(f, C.int(j)).name);
				fmt.Fprintf(w, "%s\t%d\n", gstr, i + 1);
			}
		}
	} else {
		cstr := C.CString(r.Form["q"][0]); // TODO: throw a warning where there are more
		var err C.int;
		ke := C.ke_parse(cstr, &err);
		C.free(unsafe.Pointer(cstr));
		for i := 0; i < len(bgt_files); i += 1 {
			f := bgt_files[i].f;
			for j := 0; j < int(f.n_rows); j += 1 {
				if (C.fmf_test(f, C.int(j), ke) != C.int(0)) {
					gstr := C.GoString(C.fmf_get_row(f, C.int(j)).name);
					fmt.Fprintf(w, "%s\t%d\n", gstr, i + 1);
				}
			}
		}
		C.ke_destroy(ke);
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
