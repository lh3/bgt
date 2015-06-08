#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
	next if $. == 1;
	my @t = split("\t");
	my @a = ($t[6], "altID:Z:$t[4]", "subpop:Z:$t[7]", "region:Z:$t[8]", "country:Z:$t[9]");
	push(@a, "latitude:f:$t[13]", "longtitude:f:$t[14]") if $t[13] ne '?' && $t[14] ne '?';
	print join("\t", @a), "\n";
}
