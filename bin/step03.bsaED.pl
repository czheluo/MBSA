#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$int,$dsh,$ann,$fold);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"out:s"=>\$out,
	"fold:s"=>\$fold,
	"dsh:s"=>\$dsh,
	"ann:s"=>\$ann,
			) or &USAGE;
&USAGE unless ($out);
$fold||=3;
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
my $table=ABSOLUTE_DIR("$out/01.vcf-convert");
mkdir "$out/03.bsaED" if (!-d "$out/03.bsaED");
my $int=ABSOLUTE_DIR("$out/03.bsaED");

open SH,">$dsh/03.bsaED.sh";
print SH "Rscript $Bin/bin/step03_ED.R --input $out/02.bsaRData/ --output $int && ";
print SH "perl $Bin/bin/region.pl -select $int/ED1/ED1.result -out $int/ED1/ED.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $int/ED1/ED.region -table $table/pop.vcf.table -out $int/ED1/region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $int/region.threshold.gene -i $int/ED.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $int/ED1/region.threshold.gene.kegg.stat --output $int/ED1/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $int/ED1/region.threshold.gene.go.stat --output $int/ED1/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $int/ED1/region.threshold.gene.eggnog.stat --output $int/ED1/region.threshold.gene.eggnog.stat --eggnog && ";

#print SH "Rscript $Bin/bin/step03_ED.R --input $out/02.bsaRData/ --output $int && ";
print SH "perl $Bin/bin/region.pl -select $int/ED2/ED2.result -out $int/ED2/ED.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $int/ED2/ED.region -table $table/pop.vcf.table -out $int/ED2/region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $int/region.threshold.gene -i $int/ED.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $int/ED2/region.threshold.gene.kegg.stat --output $int/ED2/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $int/ED2/region.threshold.gene.go.stat --output $int/ED2/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $int/ED2/region.threshold.gene.eggnog.stat --output $int/ED2/region.threshold.gene.eggnog.stat --eggnog && ";

#print SH "Rscript $Bin/bin/step03_ED.R --input $out/02.bsaRData/ --output $int && ";
print SH "perl $Bin/bin/region.pl -select $int/ED3/ED3.result -out $int/ED3/ED.region && ";
print SH "perl $Bin/bin/region-variant.pl -vcf $out/pop.final.vcf -select $int/ED3/ED.region -table $table/pop.vcf.table -out $int/ED3/region.variant &&";
print SH "perl $Bin/bin/region-gene.pl -a $ann -o $int/region.threshold.gene -i $int/ED.region && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $int/ED3/region.threshold.gene.kegg.stat --output $int/ED3/region.threshold.gene.kegg.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $int/ED3/region.threshold.gene.go.stat --output $int/ED3/region.threshold.gene.go.stat --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $int/ED3/region.threshold.gene.eggnog.stat --output $int/ED3/region.threshold.gene.eggnog.stat --eggnog && ";

close SH;
my $job="qsub-slurm.pl --Queue dna --Resource mem=10G --CPU 6  $dsh/03.bsaED.sh";
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
