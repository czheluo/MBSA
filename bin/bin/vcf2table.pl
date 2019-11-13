#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$P1,$P2,$B1,$B2,$Pdep,$Bdep,$popt,$Vtype,);
GetOptions(
				"help|?" =>\&USAGE,
				"vcf:s"=>\$fIn,
				"out:s"=>\$fOut,
				"Vtype:s"=>\$Vtype,
				"p1:s"=>\$P1,
				"p2:s"=>\$P2,
				"b1:s"=>\$B1,
				"b2:s"=>\$B2,
				"pdep:s"=>\$Pdep,
				"Bdep:s"=>\$Bdep,
				"popt:s"=>\$popt,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $B1 and $B2);
$Pdep||=11;
$Bdep||=10;
$popt||="F2";
$Vtype||="ALL";
open In,$fIn;
open Out,">$fOut";
my @Sample;
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	if (/^#/) {
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@info)=split(/\t/,$_);
		push @Sample,@info;
		my @out;
		push @out,join("\t","$P1.GT","$P1.AD","$P1.DP") if (defined $P1);
		push @out,join("\t","$P2.GT","$P2.AD","$P2.DP") if (defined $P2);
		push @out,join("\t","$B1.GT","$B1.AD","$B1.DP") if (defined $B1);
		push @out,join("\t","$B2.GT","$B2.AD","$B2.DP") if (defined $B2);
		print Out join("\t","Chrom\tPos\tRef\tAlt\tVtype",@out,"Ann"),"\n";
	}else{
		my (%samples,$vtype,@ann,$baseinfo,$allenum);
		$vtype=readvcf($_,\%samples,\@Sample,\@ann,\$baseinfo,\$allenum);

		if ($vtype ne $Vtype && $Vtype ne "ALL"){
			$stat{vtype}++;
			next;
		}
		if ($samples{$B1}{GT} eq "NN" || $samples{$B1}{DP} < $Bdep){
			$stat{B1}++;
			next;
		}
		if ($samples{$B2}{GT} eq "NN" || $samples{$B2}{DP} < $Bdep){
			$stat{B2}++;
			next;
		}
		if (defined $P1 && ($samples{$P1}{GT} eq "NN"||$samples{$P1}{DP} < $Pdep)){
			$stat{P1}++;
			next;
		}
		if (defined $P2 && ($samples{$P2}{GT} eq "NN"||$samples{$P2}{DP} < $Pdep)){
			$stat{P2}++;
			next;
		}
		if ($allenum > 2){
			$stat{allenum}++;
			next;
		};
		if ($popt ne "F1" && defined $P1 && defined $P2 && $samples{$P1}{GT} eq $samples{$P2}{GT}){
			$stat{PMsame}++;
			next;
		};

		my @g1=split("/",$samples{$P1}{GT});
		my @g2=split("/",$samples{$P2}{GT});
		if ($popt ne "F1" && ($g1[0] ne $g1[1]||$g2[0] ne $g2[1])){
			$stat{F2}++;
			next;
		};
		if ($popt eq "F1" && $g1[0] eq $g1[1] && $g2[0] eq $g1[1]){
			$stat{F1}++;
			next;
		}

		my @out;
		push @out,join("\t",$samples{$P1}{GT},$samples{$P1}{AD},$samples{$P1}{DP}) if (defined $P1);
		push @out,join("\t",$samples{$P2}{GT},$samples{$P2}{AD},$samples{$P2}{DP}) if (defined $P2);
		push @out,join("\t",$samples{$B1}{GT},$samples{$B1}{AD},$samples{$B1}{DP}) if (defined $B1);
		push @out,join("\t",$samples{$B2}{GT},$samples{$B2}{AD},$samples{$B2}{DP}) if (defined $B2);
		print Out join("\t",$baseinfo,$vtype,@out,join(";",@ann)),"\n";
	}
}
close In;
close Out;
#print Dumper %stat;

#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub readvcf{
	my ($line,$sample,$Sample,$anninfo,$baseinfo,$allenum)=@_;
	my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@info)=split(/\s+/,$line);
	$$baseinfo=join("\t",$chrom,$pos,$ref,$alt);
	my @alles=split(",",join(",",$ref,$alt));
	$$allenum=scalar @alles;
	if($info=~/ANN=([^\;]*)/g){
		my @ann=split(/\,/,$1);
		for (my $i=0;$i<@ann;$i++) {
			my @str=split(/\|/,$ann[$i]);
			$str[0]||="--";
			$str[1]||="--";
			$str[2]||="--";
			$str[3]||="--";
			$str[4]||="--";
			my $ann=join("|",$str[0],$str[1],$str[2],$str[3],$str[4]);
			push @{$anninfo},$ann;
		}
	}
	my %len;
	for (my $i=0;$i<@alles;$i++) {
		$len{length($alles[$i])}=1;
	}
	my $type="SNP";
	if (scalar keys %len > 1) {
		$type="INDEL";
	}
	my @format=split(/\:/,$format);
	for (my $i=0;$i<@info;$i++) {
		my @infos=split(/\:/,$info[$i]);
		for (my $j=0;$j<@infos;$j++) {
			if ($format[$j] eq "GT") {
				if ($infos[$j] =~ /\./) {
					$$sample{$$Sample[$i]}{$format[$j]}="NN";
					$$sample{$$Sample[$i]}{DP}=0;
					$$sample{$$Sample[$i]}{AD}=0;
				}else{
					my @gt=split(/\//,$infos[$j]);
					$$sample{$$Sample[$i]}{$format[$j]}=join("/",sort($alles[$gt[0]],$alles[$gt[1]]));
				}
			}
			if ($format[$j] eq "AD") {
				$$sample{$$Sample[$i]}{$format[$j]}=$infos[$j];
			}
			if ($format[$j] eq "DP") {
				$$sample{$$Sample[$i]}{$format[$j]}=$infos[$j];
			}
		}
	}
	return $type;
}
sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: czheluo\@gmail.com

Usage:
  Options:
	-vcf	<file>	input file name
	-out	<file>	output file name 
	-p1	 <str>	mut parent id name
	-p2	 <str>	wild parent id name
	-b1	<str>	mut bulk id name
	-b2	<str>	wild bulk id name
	-pdep	<str>	parent depth default 10
	-bdep	<str>	mut build depth default 10
	-popt	<str>	population type default F2
USAGE
	print $usage;
	exit;
}
