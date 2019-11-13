#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$tgl,$bs,$threshold,$dsh,$ann);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ann:s"=>\$ann,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"threshold:s"=>\$threshold,
	"tgl:s"=>\$tgl,
	"bs:s"=>\$bs,
			) or &USAGE;
&USAGE unless ($out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$bs||=30;
$tgl||=2000;
$threshold||=0.999;
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
my $table=ABSOLUTE_DIR("$out/01.vcf-convert");
mkdir "$out/06.bsaBSR" if (!-d "$out/06.bsaBSR");
my $out3=ABSOLUTE_DIR("$out/06.bsaBSR");

open SH,">$dsh/06.bsaBSR.sh";
print SH "Rscript $Bin/bin/step06_BSR.R --input $out/02.bsaRData --output $out3 --tgl $tgl --bs $bs --threshold $threshold && ";
print SH "perl $Bin/bin/region.pl -select $out3/BSR.result -out $out3/BSR.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $out3/BSR.region -table $table/pop.vcf.table -out $out3/BSR.region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out3/region.threshold.gene -i $out3/BSR.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/region.threshold.gene.kegg.stat --output $out3/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/region.threshold.gene.go.stat --output $out3/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/region.threshold.gene.eggnog.stat --output $out3/region.threshold.gene.eggnog.stat --eggnog && ";
close SH;
my $job="qsub-slurm.pl --Queue dna --Resource mem=5G --CPU 1  $dsh/06.bsaBSR.sh";
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
Contact:        czheluo\@gmail.com
Script:			$Script
Description:
	
	eg:
	perl  

Usage:
  Options:
	-int input dir
	-out output dir 
	-fold 
	-ann annotation summary
	-dsh work shell 
    -h      Help

USAGE
        print $usage;
        exit;
}
