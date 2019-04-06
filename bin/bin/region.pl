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
my $start = -1;
my $CHR = "chr1";
my ($POS,%region,$THRE);
my $PVAL = 0;
my $END = 0;
while (<IN>) {
	chomp;
	#s/"//g;
	next if ($_ eq "" || /^$/ || /^Chr/);
	my($chr,$pos,$pval,$thre)=split/\s+/;
	my $end=-2;
	if ($CHR ne $chr) {
		if ($PVAL > $thre) {
			$end = $POS;
			$region{$CHR}{join("\t",$start,$end)}=1;
		}
		if ($pval > $thre) {
			$start = $pos;
			$END = -1;
		}else{
			$start = -1;
		}
	}else{
		if ($start eq "-1") {
			if ($pval > $thre) {
				$start = $pos;
				$END = -1;
			}
		}else{
			if ($END eq "-1") {
				if ($pval < $thre) {
					$end = $POS;
					$region{$chr}{join("\t",$start,$end)}=1;
					$start = -1;
				}
			}
		}
	}

	$CHR = $chr;
	$POS = $pos;
	$PVAL = $pval;
	$THRE = $thre;
	#print "$start\t$end\t666\t$END\n";
}
close IN;
#print "$CHR\t$start\t$POS\t$PVAL\t$THRE\t$END";
#$region{$CHR}{join("\t",$start,$POS)}=1 if ($PVAL >= $THRE && $END eq "-1");

#print Dumper \%region;die;
open OUT,">",$out;
print OUT "#chr\tstart\tend\n";
foreach my $chr(sort keys %region) {
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($start,$end)=split(/\s+/,$region);
		next if ($start eq "-1" || $end eq "1");

		print OUT "$chr\t$region\n";
	}
}
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        minghao.zhang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	-select the methods result file name
	-out output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
