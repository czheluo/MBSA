#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin1,$fin2,$fout,$dp1,$dp2,,$min);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"int1:s"=>\$fin1,	
	"int2:s"=>\$fin2,
	"wdp:s"=>\$dp1,
	"mdp:s"=>\$dp2,
	"out:s"=>\$fout,
			) or &USAGE;
&USAGE unless ($fout);

open In,$fin2;
open Out,">$fout/pop.new.table";
while (<In>) {
	chomp;
	if ($_ =~ "CHROM") {
		print Out "$_\n";
	}else{
		my($all)=split/\s+/,$_,1;
		my @alls=split/\s+/,$all;
		#print $alls[-5];;die;
		my $wdp=$alls[-5];
		my $mdp=$alls[-2];
		next if ($wdp < $dp1);
		next if ($mdp < $dp2);
		print Out "$_\n";
	}
	
}
close In;
close Out;
open IN,$fin1;
open OUT,">$fout/pop.new.index";
while (<IN>) {
	chomp;
	if ($_ =~ "#chr") {
		print OUT "$_\n";
	}else{
		my($all)=split/\s+/,$_,1;
		my @alls=split/\s+/,$all;
		#print $alls[-5];;die;
		my $wdp=$alls[-4]+$alls[-2];
		my $mdp=$alls[-1]+$alls[-3];
		next if ($wdp < $dp1);
		next if ($mdp < $dp2);
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

	eg: perl $Script -int2 pop.tables -int1 pop.bsa.index -out ./ -wdp 70 -mdp 70
	

Usage:
  Options:
	"int1:s"=>\$fin1,	
	"int2:s"=>\$fin2,
	"wdp:s"=>\$dp1,
	"mdp:s"=>\$dp2,
	"out:s"=>\$fout,
	-h         Help

USAGE
        print $usage;
        exit;
}
