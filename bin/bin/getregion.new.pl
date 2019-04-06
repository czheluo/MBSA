#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($select,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"select:s"=>\$select,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($out);
open IN,$select;
open OUT,">",$out;
print OUT "#chr\tstart\tend\n";
my $start = -1;
my $end = -1;
my $CHR = "chr1";
my $POS;
while (<IN>) {
	chomp;
	#s/"//g;
	next if ($_ eq "" || /^$/ || /Chr/);
	my($chr,$pos,$pval,$thre)=split/\s+/;
	if ($CHR ne $chr) {
		if ($pval > $thre) {
			$start = $pos;
		}else{
			$start = -1;
		}
	}else{
		if ($start eq "-1") {
			if ($pval > $thre) {
				$start = $pos;
				$end = -1;
			}
		}else{
			if ($end eq "-1") {
				if ($pval < $thre) {
					$end = $POS;
					print OUT "$CHR\t$start\t$end\n";
					$start = -1;
				}else{
					$end = $pos;
					print OUT "$chr\t$start\t$end\n";
				}
			}
		}
	}
	$CHR = $chr;
	$POS = $pos;
}
close IN;
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"select:s"=>\$select,
	"out:s"=>\$out,
  -h         Help

USAGE
        print $usage;
        exit;
}
