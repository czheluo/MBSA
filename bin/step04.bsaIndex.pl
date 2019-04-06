#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$ann,$pta,$ptb,$bs,$ws,$popt,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ann:s"=>\$ann,
	"out:s"=>\$out,
	"pta:s"=>\$pta,
	"ptb:s"=>\$ptb,
	"bs:s"=>\$bs,
	"ws:s"=>\$ws,
	"popt:s"=>\$popt,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$popt||="F2";
$pta||=95;
$ptb||=99;
$bs||=30;
$ws||=1e6;
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
my $table=ABSOLUTE_DIR("$out/01.vcf-convert");
mkdir "$out/04.bsaIndex" if (!-d "$out/04.bsaIndex");
my $out3=ABSOLUTE_DIR("$out/04.bsaIndex");

open SH,">$dsh/04.bsaindex.sh";
print SH "Rscript $Bin/bin/step04_fitIndex.R --input $out/02.bsaRData --output $out3  --pta $pta --ptb $ptb --ws $ws --bs $bs && ";
print SH "perl $Bin/bin/region.pl -select $out3/Index1/Index1.result -out $out3/Index1/Index.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $out3/Index1/Index.region -table $table/pop.vcf.table -out $out3/Index1/Index.region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out3/Index1/region.threshold.gene -i $out3/Index1/Index.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index1/region.threshold.gene.kegg.stat --output $out3/Index1/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index1/region.threshold.gene.go.stat --output $out3/Index1/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index1/region.threshold.gene.eggnog.stat --output $out3/Index1/region.threshold.gene.eggnog.stat --eggnog && ";

print SH "perl $Bin/bin/region.pl -select $out3/Index2/Index2.result -out $out3/Index2/Index.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $out3/Index2/Index.region -table $table/pop.vcf.table -out $out3/Index2/Index.region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out3/Index2/region.threshold.gene -i $out3/Index2/Index.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index2/region.threshold.gene.kegg.stat --output $out3/Index2/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index2/region.threshold.gene.go.stat --output $out3/Index2/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index2/region.threshold.gene.eggnog.stat --output $out3/Index2/region.threshold.gene.eggnog.stat --eggnog && ";

print SH "perl $Bin/bin/region.pl -select $out3/Index3/Index3.result -out $out3/Index3/Index.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $out3/Index3/Index.region -table $table/pop.vcf.table -out $out3/Index3/Index.region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out3/Index3/region.threshold.gene -i $out3/Index3/Index.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index3/region.threshold.gene.kegg.stat --output $out3/Index3/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index3/region.threshold.gene.go.stat --output $out3/Index3/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index3/region.threshold.gene.eggnog.stat --output $out3/Index3/region.threshold.gene.eggnog.stat --eggnog && ";

print SH "perl $Bin/bin/region.pl -select $out3/Index4/Index4.result -out $out3/Index4/Index.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $out3/Index4/Index.region -table $table/pop.vcf.table -out $out3/Index4/Index.region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out3/Index4/region.threshold.gene -i $out3/Index4/Index.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index4/region.threshold.gene.kegg.stat --output $out3/Index4/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index4/region.threshold.gene.go.stat --output $out3/Index4/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Index4/region.threshold.gene.eggnog.stat --output $out3/Index4/region.threshold.gene.eggnog.stat --eggnog && ";

close SH;
my $job="qsub-slurm.pl --Queue dna --Resource mem=3G --CPU 1  $dsh/04.bsaindex.sh";
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
	-ann pop.summary or anno.summary
	-out out dir 
	-bs    bulk size default was 30
    -ws window size default was 1M
	-pta pemutation test confidence interval
	-ptb pemutation test confidence interval

    -bs Set N to be the number of individuals in the mutant pool

    -h      Help

USAGE
        print $usage;
        exit;
}
