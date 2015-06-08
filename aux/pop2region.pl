#!/usr/bin/env perl

use strict;
use warnings;

my %h = (CHB=>'EastAsia', JPT=>'EastAsia', CHS=>'EastAsia', CDX=>'EastAsia', KHV=>'EastAsia', CHD=>'EastAsia',
		 CEU=>'WestEurasia', TSI=>'WestEurasia', GBR=>'WestEurasia', FIN=>'WestEurasia', IBS=>'WestEurasia',
		 YRI=>'Africa', LWK=>'Africa', GWD=>'Africa', MSL=>'Africa', ESN=>'Africa',
		 ASW=>'America', ACB=>'America', MXL=>'America', PUR=>'America', CLM=>'America', PEL=>'America',
		 GIH=>'SouthAsia', PJL=>'SouthAsia', BEB=>'SouthAsia', STU=>'SouthAsia', ITU=>'SouthAsia');

while (<>) {
	chomp;
	if (/population:Z:(\S+)/ && defined($h{$1})) {
		print "$_\tregion:Z:$h{$1}\tsource:Z:1000G\n";
	}
}
