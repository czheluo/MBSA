#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
my $resultdir = "$out/result";
$out = "$out/work_flow";
`cp $out/variant.region.stat $resultdir/Table/3-12.xls`;
`perl $Bin/gatherdata.pl -i $out/region.threshold.gene.total -o $resultdir/Table`;
`cut -d "," -f1,3-15 $out/qtl.csv >$resultdir/Result/qtl.region.xls`;
`grep "\@chr" $out/region.threshold.variant.total|sed 's/@//g' > $resultdir/Table/3-13.xls `;
`cp $out/region*.total $resultdir/Result`;
`cp $out/region*.eff $resultdir/Result`;
`cp $out/*.detail $resultdir/Result`;
`cp $out/*.pdf $resultdir/Figure`;
`cp $out/*.png $resultdir/Figure`;
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"help|?" =>\&USAGE,
	"out:s"=>\$out, output dir
USAGE
        print $usage;
        exit;
}
