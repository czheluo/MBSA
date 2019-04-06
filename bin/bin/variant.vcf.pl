#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$vcf);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"G:s"=>\$fin,
	"vcf:s"=>\$vcf,
	"o:s"=>\$fout,
			) or &USAGE;
&USAGE unless ($fout);
open IN,$fin;
my %hash;
while (<IN>) {
	chomp;
	s/\"//g;
	next if ($_ eq "" || /^$/ || /CHROM/);
	my ($chr,$pos,undef)=split(/\t/,$_);
	my $id = $chr."_".$pos;
	$hash{$id}=1;
}
close IN;
open IN,$vcf;
open OUT,">$fout/variant.region.vcf";
my %stat;
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/ || /##/);
	if (/#/) {
		print OUT "$_\n";
	}else{
		my ($chr,$pos,$id,$ref,$alt,undef)=split(/\t/,$_);
		if (exists $hash{$id}) {
			print OUT "$_\n";
			my @alle=split(",",join(",",$ref,$alt));
			my %len;
			for (my $i=0;$i<@alle;$i++) {
				$len{length($alle[$i])}=1;
			}
			if (scalar keys %len>1) {
				$stat{$chr}{indel}++;
			}else{
				$stat{$chr}{snp}++;
			}
		}
	}
}
close IN;
close OUT;
open OUT ,">$fout/variant.region.stat";
print OUT "chr\tSNP number\tIndel number\n";
foreach my $chr (sort keys %stat) {
	$stat{$chr}{snp}||=0;
	$stat{$chr}{indel}||=0;
	print OUT "$chr\t$stat{$chr}{snp}\t$stat{$chr}{indel}\n";
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
