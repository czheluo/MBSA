#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$ann,$bs,$ws,$popt,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ann:s"=>\$ann,
	"out:s"=>\$out,
	"ws:s"=>\$ws,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$bs||=30;
$ws||=1e6;
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
my $table=ABSOLUTE_DIR("$out/01.vcf-convert");
mkdir "$out/05.bsaGprime" if (!-d "$out/05.bsaGprime");
my $out3=ABSOLUTE_DIR("$out/05.bsaGprime");

open SH,">$dsh/05.bsaGprime.sh";
print SH "Rscript $Bin/bin/step05_Gprime.R --input $out/02.bsaRData --output $out3  --ws $ws --bs $bs && ";
print SH "perl $Bin/bin/region.pl -select $out3/Gprime1/Gprime1.result -out $out3/Gprime1/Gprime.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $out3/Gprime1/Gprime.region -table $table/pop.vcf.table -out $out3/Gprime1/Gprime.region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out3/Gprime1/region.threshold.gene -i $out3/Gprime1/Gprime.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime1/region.threshold.gene.kegg.stat --output $out3/Gprime1/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime1/region.threshold.gene.go.stat --output $out3/Gprime1/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime1/region.threshold.gene.eggnog.stat --output $out3/Gprime1/region.threshold.gene.eggnog.stat --eggnog && ";

print SH "perl $Bin/bin/region.pl -select $out3/Gprime2/Gprime2.result -out $out3/Gprime2/Gprime.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $out3/Gprime2/Gprime.region -table $table/pop.vcf.table -out $out3/Gprime2/Gprime.region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out3/Gprime2/region.threshold.gene -i $out3/Gprime2/Gprime.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime2/region.threshold.gene.kegg.stat --output $out3/Gprime2/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime2/region.threshold.gene.go.stat --output $out3/Gprime2/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime2/region.threshold.gene.eggnog.stat --output $out3/Gprime2/region.threshold.gene.eggnog.stat --eggnog && ";

print SH "perl $Bin/bin/region.pl -select $out3/Gprime3/Gprime3.result -out $out3/Gprime3/Gprime.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $out3/Gprime3/Gprime.region -table $table/pop.vcf.table -out $out3/Gprime3/Gprime.region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out3/Gprime3/region.threshold.gene -i $out3/Gprime3/Gprime.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime3/region.threshold.gene.kegg.stat --output $out3/Gprime3/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime3/region.threshold.gene.go.stat --output $out3/Gprime3/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime3/region.threshold.gene.eggnog.stat --output $out3/Gprime3/region.threshold.gene.eggnog.stat --eggnog && ";

print SH "perl $Bin/bin/region.pl -select $out3/Gprime4/Gprime4.result -out $out3/Gprime4/Gprime.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $out3/Gprime4/Gprime.region -table $table/pop.vcf.table -out $out3/Gprime4/Gprime.region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $out3/Gprime4/region.threshold.gene -i $out3/Gprime4/Gprime.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime4/region.threshold.gene.kegg.stat --output $out3/Gprime4/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime4/region.threshold.gene.go.stat --output $out3/Gprime4/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out3/Gprime4/region.threshold.gene.eggnog.stat --output $out3/Gprime4/region.threshold.gene.eggnog.stat --eggnog && ";

close SH;
my $job="qsub-slurm.pl --Queue dna --Resource mem=3G --CPU 1  $dsh/05.bsaGprime.sh";
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
	-ann pop.summary or anno.summary
	-out out dir 
	-bs    bulk size default was 30
    -ws window size default was 1M
    -h      Help

USAGE
        print $usage;
        exit;
}
