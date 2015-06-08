#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
	next if /^#/;
	chomp;
	my @t = split("\t");
	my ($chr, $st, $en);
	if ($t[1] =~ /^(\S+):(\d+)(-(\d+))?/) {
		$chr = $1;
		$st = int($2);
		$en = defined($3)? int($4) : $st;
	}
	my $al;
	my $is_sym = 0;
	if ($t[2] eq '-') { # deletion
		$al = "$chr:$st:" . ($en - $st + 1) . ":";
	} elsif ($st == $en && length($t[2]) == 1) { # SNP
		$al = "$chr:$st:1:$t[2]";
	} elsif ($en - $st == 1) { # insertion
		$al = "$chr:$en:0:$t[2]";
	} else {
		$al = "$chr:$st:" . ($en - $st + 1) . ":<$t[2]>";
		$is_sym = 1;
	}
	next if $is_sym;

	my @s = ($al);
	my @u = split(",", $t[6]);
	push(@s, "rsID:Z:$t[12]") if $t[12] ne '-';
	for my $x (@u) {
		push(@s, "effect:Z:$x");
	}
	if ($t[13] =~ /IMPACT=([^\s;]+)/) {
		push(@s, "impact:Z:$1");
	}
	if ($t[13] =~ /SYMBOL=([^\s;]+);SYMBOL_SOURCE=HGNC;.*BIOTYPE=([^\s;]+)/) {
		push(@s, "gene:Z:$1", "biotype:Z:$2");
	}
	if ($t[4] ne '-' && $t[5] ne '-') {
		push(@s, "featureID:Z:$t[4]", "featureType:Z:$t[5]");
	}
	push(@s, "CCDS:Z:$1") if $t[13] =~ /;CCDS=([^\s;]+)/;
	push(@s, "CDSpos:i:$t[8]") if $t[8] ne '-';
	if ($t[13] =~ /DISTANCE=(\d+);STRAND=(\d+)/) {
		push(@s, "distance:i:$1", "strand:i:$2");
	}
	push(@s, "codon:Z:$t[11]") if $t[10] ne '-';
	if ($t[13] =~ /SIFT=([^\s;()]+)\(([\d.]+)\)/) {
		push(@s, "SIFT:Z:$1");
		push(@s, "SIFT_p:f:$2") if defined($2);
	}
	if ($t[13] =~ /PolyPhen=([^\s;()]+)\(([\d.]+)\)/) {
		push(@s, "PolyPhen:Z:$1");
		push(@s, "PolyPhen_p:f:$2") if defined($2);
	}
	print(join("\t", @s), "\n");
}
