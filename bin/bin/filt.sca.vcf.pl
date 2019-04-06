#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($table,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"table:s"=>\$table,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($out);
open IN,$table;
open OUT,">$out";
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/ );
	if (/^Chrom/) {
		print OUT "$_\n";
	}else{
		my ($chr,undef)=split(/\t/,$_);
		next if ($chr =~ /sca/);
		print OUT "$_\n";
	}
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
	eg:
	perl filt.sca.vcf.pl -table pop.final.vcf.table -out pop.vcf.table

Usage:
  Options:
  -table	<file>	input table file name
  -out	<file>	out no sca table file name 
  -h         Help

USAGE
        print $usage;
        exit;
}
