#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$ann,$pid,$bid,$step,$stop,$pta,$ptb,$bs,$ws,$tgl,$popt);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"ann:s"=>\$ann,
	"pid:s"=>\$pid,
    	"bid:s"=>\$bid,
	"out:s"=>\$out,
	"pta:s"=>\$pta,
	"ptb:s"=>\$ptb,
	"bs:s"=>\$bs,
	"ws:s"=>\$ws,
	"popt:s"=>\$popt,
	"tgl:s"=>\$tgl,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($vcf and $out);
$popt||="F2";
$pta||=95;
$ptb||=99;
$bs||=30;
$ws||=1e6;
$tgl||=2000;
$step||=1;
$stop||=0;
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf); 
my $dsh="$out/work_sh";
mkdir $dsh if (!-d $dsh);
open Log,">$out/work_sh/BSA.$BEGIN_TIME.log";
if ($step == 1) {
	print Log "########################################\n";
	print Log "vcf-convert \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.vcf-convert.pl -vcf $vcf  -out $out -pid $pid -bid $bid -dsh $dsh\n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "vcf-convert Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 2) {
	print Log "########################################\n";
	print Log "BSA-DATA \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step02.bsaRData.pl -out $out -pid $pid -bid $bid -ws $ws -dsh $dsh\n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "prepare data Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 3) {
	print Log "########################################\n";
	print Log "bsa-ED \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step03.bsaED.pl -ann $ann -out $out -dsh $dsh \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Calculated ED Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 4) {
	print Log "########################################\n";
	print Log "BSA-Index \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step04.bsaIndex.pl -ann $ann -popt $popt -out $out -pta $pta -ptb $ptb -bs $bs -ws $ws -dsh $dsh \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Calculated Index  Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 5 ) {
	print Log "########################################\n";
	print Log "BSA-Gprime\n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step05.bsaGprime.pl -ann $ann -ws $ws -out $out -dsh $dsh \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Calculated G¡¯ Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 6) {
	print Log "########################################\n";
	print Log "BSA-BSR \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step06.bsaBSR.pl -ann $ann -tgl $tgl -bs $bs -out $out -dsh $dsh \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "calculated BSR Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
close Log;
#`cp $out/variant.region.stat $resultdir/Table/3-12.xls`;
#`perl $Bin/gatherdata.pl -i $out/region.threshold.gene.total -o $resultdir/Table`;
#`cut -d "," -f1,3-15 $out/qtl.csv >$resultdir/Result/qtl.region.xls`;
#`grep "\@chr" $out/region.threshold.variant.total|sed 's/@//g' > $resultdir/Table/3-13.xls `;
##`cp $out/region*.total $resultdir/Result`;
#`cp $out/region*.eff $resultdir/Result`;
#`cp $out/*.stat $resultdir/Result`;
#`cp $out/*.pdf $resultdir/Figure`;
#`cp $out/*.png $resultdir/Figure`;

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
Contact:       czheluo\@gmail.com
Script:			$Script
Description: 
	Multple methods for BSA Pipeline 
Usage:
  Options:
	-vcf   pop.final.vcf
	-ann pop.summary or anno.summary
	-ref ref.fa file or if you want to use the GATK to generate TABLE (default is NO)
	-pid  parental name A,B
	-bid  bulk name C,D
	-out out dir 
	-bs    bulk size default was 30
    	-ws window size default was 1M
	-pta pemutation test confidence interval
	-ptb pemutation test confidence interval
	-tgl  total genetic length (default was 2000)
    	-bs Set N to be the number of individuals in the mutant pool
	-step which step you want 
	-stop control the step 
    	-h      Help

USAGE
        print $usage;
        exit;
}
