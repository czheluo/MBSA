#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$pid,$bid,$popt,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"pid:s"=>\$pid,
    "bid:s"=>\$bid,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($vcf and $out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf); 
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
mkdir "$out/01.vcf-convert" if (!-d "$out/01.vcf-convert");
my $out=ABSOLUTE_DIR("$out/01.vcf-convert");

my @bid=split(/\,/,$bid);
my @pid=split(/\,/,$pid);
open SH,">$dsh/01.vcf.convert.sh";
if (scalar @bid == 1) {
	print SH "perl $Bin/bin/vcf2table.pl -vcf $vcf -B1 $bid -out $out/pop.final.vcf.table ";
	if (scalar @pid==1) {
		print SH "-P1 $pid &&";
	}elsif(scalar @pid==2)  {
		print SH "-P1 $pid[0] -P2 $pid[1] &&";
	}else{
		print SH "&&";
	}
}else{
	print SH "perl $Bin/bin/vcf2table.pl -vcf $vcf  -B1 $bid[0] -B2 $bid[1] -out $out/pop.final.vcf.table ";
	if (scalar @pid==1) {
		print SH "-P1 $pid &&";
	}elsif(scalar @pid==2)  {
		print SH "-P1 $pid[0] -P2 $pid[1] &&";
	}else{
		print SH "&&";
	}
}
print SH "perl $Bin/bin/filt.sca.vcf.pl -table $out/pop.final.vcf.table -out $out/pop.vcf.table && rm $out/pop.final.vcf.table ";
#print SH "java -Xmx50G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T VariantsToTable -R $ref -V $vcf -F CHROM -F POS -F REF -F ALT -GF AD -GF DP -GF GQ -GF PL -o $out/pop.vcf.table";
close SH;
my $job="qsub-slurm.pl --Queue dna --Resource mem=4G --CPU 2  $dsh/01.vcf.convert.sh";
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
	-vcf   pop.final.vcf
	-pid  parental name A,B
	-bid  bulk name C,D
	-out out dir 
	-dsh work shell
    -h      Help

USAGE
        print $usage;
        exit;
}
