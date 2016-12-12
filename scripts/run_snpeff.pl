#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use File::Path qw(make_path remove_tree);
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;


my $jar_snpeff = $Setting->{jar_snpeff} ;
my $dir_snpeff = $Setting->{dir_snpeff} ;
my $genome_build = $Setting->{genome_snpeff} ;

# Options handling
my $opts = {} ;
getopts("i:o:l:j:n:d:p:", $opts) ;

my $jobid = $opts->{j} ;
my $threads_active = $opts->{n} || 16 ;
my $file_vcfgz_input  = $opts->{i} ;
my $file_vcfgz_output = $opts->{o} ;
my $file_geminidb = $opts->{d} ;
my $file_ped = $opts->{p}  ;
my $dir_log = $opts->{l} ;
my $file_log = "$dir_log/running_snpeff.$jobid.log" ;

my $file_vcf_snpeff_tmp = "$dir_log/snpeff_tmp.vcf" ;
my $file_db_snpeff_tmp = "$dir_log/snpeff_tmp.db" ;
my $file_gemini_dump  = "$dir_log/snpeff_geminidump.txt" ;
my $file_vcfgz_snpeff_dump_fromdb = "$dir_log/snpeff_gdump.vcf.gz" ;


# Step 1. snpeff
open (my $LOG, ">$file_log") ;
my $cmd_snpeff =  "java -jar $jar_snpeff -dataDir $dir_snpeff $genome_build $file_vcfgz_input > $file_vcf_snpeff_tmp" ;
logNrun($cmd_snpeff) ;

bgzip($file_vcf_snpeff_tmp, "$file_vcf_snpeff_tmp.gz") ;
tabix_vcf("$file_vcf_snpeff_tmp.gz") ;

# Step 2. gemini db
my $cmd_gemini_load = "gemini load " ;
$cmd_gemini_load .= " --cores $threads_active " ;
$cmd_gemini_load .= " -t snpEff " ;
$cmd_gemini_load .= " -v $file_vcf_snpeff_tmp.gz " ;
#$cmd_gemini_load .= " -p $file_ped " ;
$cmd_gemini_load .= " $file_db_snpeff_tmp 2> $dir_log/stderr.geminiload.snpeff.log" ;

logNrun($cmd_gemini_load) ;

# Step 3. gquery dump
my $cmd_gquery_dump = "gemini query --header -q 'select chrom,start,end,* from variants' $file_db_snpeff_tmp> $file_gemini_dump " ;
logNrun($cmd_gquery_dump) ;
bgzip($file_gemini_dump , "$file_gemini_dump.gz" ) ;
my $cmd_tabix = "tabix -s 1 -b 2 -e 3 -S 1 -f $file_gemini_dump.gz" ;
logNrun($cmd_tabix) ;

my $header_dump = `head -1 $file_gemini_dump ` ;
chomp($header_dump) ;
my @header_col = split(/\t/,$header_dump) ;

#skip first 4
shift @header_col ;
shift @header_col ;
shift @header_col ;

#bgzip($file_gemini_dump , "$file_gemini_dump.gz") ;
#tabix_vcf($file_gemini_dump) ;

# Step 4. vcf-annotate LOF NMD
my $cmd_vannot_gdump = "cat $file_vcf_snpeff_tmp | vcf-annotate -a $file_gemini_dump.gz " ;
$cmd_vannot_gdump .= " -c " . join("," , (qw"CHROM FROM TO" , map {"INFO/snpeff_$_" } @header_col) ) ; #zzzz

foreach my $col (@header_col) {
    $cmd_vannot_gdump .= " -d key=INFO,ID=snpeff_$col,Number=1,Type=String,Description='$col by snpEff'" ;
}

$cmd_vannot_gdump .= " 2> $dir_log/stderr_vannot_snpeff_dump.log" ;
$cmd_vannot_gdump .= " | bgzip -@ $threads_active -c > $file_vcfgz_snpeff_dump_fromdb" ;

logNrun($cmd_vannot_gdump) ;
tabix_vcf($file_vcfgz_snpeff_dump_fromdb) ;

# Step 5. snp-eff
# Step 6. gannot

my $cmd_gannot = "gemini annotate -f $file_vcfgz_snpeff_dump_fromdb -a extract " ;
$cmd_gannot .= " -e " . join("," , map {"snpeff_$_"} @header_col) . " " ;
$cmd_gannot .= " -t " . join("," , map {"text"} @header_col) . " " ;
$cmd_gannot .= " -c " . join("," , map {"snpeff_$_"} @header_col) . " " ;
$cmd_gannot .= " -o " . join("," , map {"first"} @header_col) . " " ;
$cmd_gannot .= " $file_geminidb 2> $dir_log/stderr.gannot_snpeff.log" ;

logNrun($cmd_gannot) ;


sub tabix_vcf {
    my $file = shift ;

    if ($file !~ /gz$/) {
        my $file_orig = $file ;
        $file = "$file_orig.gz" ;
        bgzip($file_orig,$file) ;
    }

    my $cmd_tabix_vcf = "tabix -p vcf -f $file " ;
    logNrun($cmd_tabix_vcf) ;
}

sub bgzip {
    my $file_src = shift ;
    my $file_tgt = shift ;

    my $cmd_bgzip = "bgzip -@ $threads_active -c $file_src > $file_tgt" ;
    logNrun($cmd_bgzip) ;
}

sub timelog {
    my $text = shift ;
    my $time = localtime(time) ;

    return "$time\t$text" ;
}

sub logit {
    my $msg = shift ;
    my $fh = shift ;

    print $fh "$msg\n" ;
}

sub tlog {  # timelog + logit
    my $msg = shift ;
    logit(timelog("$msg"),$LOG) ;
}

sub logNrun {
    my $cmd = shift ;

    tlog($cmd) ;
    `$cmd` ;
}

sub getts {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    return sprintf ( "%04d%02d%02d%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
}



