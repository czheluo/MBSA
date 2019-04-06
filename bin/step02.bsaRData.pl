#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$pid,$bid,$ws,$popt,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"pid:s"=>\$pid,
    "bid:s"=>\$bid,
	"out:s"=>\$out,
	"ws:s"=>\$ws,
	"popt:s"=>\$popt,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$popt||="F2";
$ws||=1e6;
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
my @bid=split(/\,/,$bid);
my @pid=split(/\,/,$pid);
my $table=ABSOLUTE_DIR("$out/01.vcf-convert");
mkdir "$out/02.bsaRData" if (!-d "$out/02.bsaRData");
my $int=ABSOLUTE_DIR("$out/02.bsaRData");

open SH,">$dsh/02.bsaRData.sh";
if (scalar @bid == 1) {
	print SH "Rscript $Bin/bin/step01_Index.R --input $table/pop.vcf.table --output $int --hb $bid --pop $popt --ws $ws ";
	if (scalar @pid==1) {
		print SH "--wid $pid &&";
	}elsif(scalar @pid==2) {
		print SH "--wid $pid[0] --mid $pid[1] &&";
	}else{
		print SH "&&";
	}
}else{
	print SH "Rscript $Bin/bin/step01_Index.R --input $table/pop.vcf.table --output $int --hb $bid[0] --lb $bid[1] --pop $popt --ws $ws ";
	if (scalar @pid==1) {
		print SH "--wid $pid &&";
	}elsif(scalar @pid==2) {
		print SH "--wid $pid[0] --mid $pid[1] &&";
	}else{
		print SH "&&";
	}
}
close SH;
my $job="qsub-slurm.pl --Resource mem=4G --CPU 3  $dsh/02.bsaRData.sh";
`$job`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
        my $cur_dir=`pwd`;chomp($cur_dir);
        my ($in)=@_;
        my $return="";
        if(-f $in){
                my $dir=dirname($in);
                my $file=basename($in);
                chdir $dir;$dir=`pwd`;chomp $dir;
                $return="$dir/$file";
        }elsif(-d $in){
                chdir $in;$return=`pwd`;chomp $return;
        }else{
                warn "Warning just for file and dir \n$in";
                exit;
        }
        chdir $cur_dir;
        return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl  

Usage:
  Options:
	-table   pop.vcf.table
	-pid  parental name A,B
	-bid  bulk name C,D
	-out out dir 
    -ws window size default was 1M
	-popt the type of population 
    -h      Help

USAGE
        print $usage;
        exit;
}
