#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use File::Path qw(make_path remove_tree);
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use Storable ;

my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

my $jar_snpeff = $Setting->{jar_snpeff} ;
my $dir_snpeff = $Setting->{dir_snpeff} ;
my $genome_build = $Setting->{genome_snpeff} ;

# Options handling
my $opts = {} ;
getopts("i:o:j:n:d:", $opts) ;

my $jobid             =	$opts->{j} ;
my $file_vcfgz_input  = $opts->{i} ;
my $file_vcfgz_output = $opts->{o} ;
my $file_geminidb     = $opts->{d} ;
my $threads_active    =	$opts->{n} || 16 ;

my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

my $file_log = "$dir_log/running_snpeff_$jobid.log" ;

my $file_vcf_snpeff_tmp = "$dir_work/snpeff_tmp.vcf" ;
my $file_db_snpeff_tmp = "$dir_work/snpeff_tmp.db" ;
my $file_gemini_dump  = "$dir_work/snpeff_geminidump.txt" ;
my $file_gemini_dump_tmp  = "$dir_work/snpeff_geminidump_tmp.txt" ;
#my $file_vcfgz_snpeff_dump_fromdb = "$dir_results_output/snpeff_gdump.vcf.gz" ;

my $file_hashref_vannot_tmp = "$dir_work/vannot_tmp.hashref" ;

my $log = VarSelect::Log->new(file => $file_log) ;

# Step 1. snpeff
open (my $LOG, ">$file_log") ;
my $cmd_snpeff =  "java -jar $jar_snpeff -dataDir $dir_snpeff $genome_build $file_vcfgz_input > $file_vcf_snpeff_tmp 2> $dir_log/stderr_snpEff_jar_$jobid.log" ;
$log->andRun($cmd_snpeff) ;

bgzip($file_vcf_snpeff_tmp, "$file_vcf_snpeff_tmp.gz") ;
tabix_vcf("$file_vcf_snpeff_tmp.gz") ;

# Step 2. gemini db
my $cmd_gemini_load = "gemini load " ;
$cmd_gemini_load .= " --cores $threads_active " ;
$cmd_gemini_load .= " -t snpEff " ;
$cmd_gemini_load .= " -v $file_vcf_snpeff_tmp.gz " ;
$cmd_gemini_load .= " $file_db_snpeff_tmp 2> $dir_log/stderr.geminiload.snpeff.log" ;

$log->andRun($cmd_gemini_load) ;

# Step 3. gquery dump
my $cmd_gquery_dump = "gemini query --header -q 'select chrom,start+1,REF,ALT,* from variants' $file_db_snpeff_tmp> $file_gemini_dump 2> $dir_log/stderr_snpeff_gdump_$jobid.log " ;
$log->andRun($cmd_gquery_dump) ;

# dirty hack for chromosome chrXX
my $cmd_dirtyhack_chr = "sed -i 's/^chr//' $file_gemini_dump > $file_gemini_dump_tmp; cat $file_gemini_dump_tmp >> $file_gemini_dump" ;
$log->andRun($cmd_dirtyhack_chr) ;

bgzip($file_gemini_dump , "$file_gemini_dump.gz" ) ;
my $cmd_tabix = "tabix -s 1 -b 2 -e 2 -S 1 -f $file_gemini_dump.gz" ;
$log->andRun($cmd_tabix) ;

my $header_dump = `head -1 $file_gemini_dump ` ;
chomp($header_dump) ;
my @header_col = split(/\t/,$header_dump) ;

#skip first 3 cols (chrom start end)
shift @header_col ;
shift @header_col ;
shift @header_col ;
shift @header_col ;

#bgzip($file_gemini_dump , "$file_gemini_dump.gz") ;
#tabix_vcf($file_gemini_dump) ;

# Step 4. vcf-annotate
my $cmd_vannot_gdump = "cat $file_vcf_snpeff_tmp | vcf-annotate -a $file_gemini_dump.gz " ;
$cmd_vannot_gdump .= " -c " . join("," , (qw"CHROM POS REF ALT" , map {"INFO/snpeff_$_" } @header_col) ) ; #zzzz

foreach my $col (@header_col) {
    $cmd_vannot_gdump .= " -d key=INFO,ID=snpeff_$col,Number=1,Type=String,Description='$col by snpEff'" ;
}

$cmd_vannot_gdump .= " 2> $dir_log/stderr_vannot_snpeff_dump_$jobid.log" ;
$cmd_vannot_gdump .= " | bgzip -@ $threads_active -c > $file_vcfgz_output" ;

$log->andRun($cmd_vannot_gdump) ;
tabix_vcf($file_vcfgz_output) ;

## Step 5. gannot
#
#my $cmd_gannot = "gemini annotate -f $file_vcfgz_output -a extract " ;
#$cmd_gannot .= " -e " . join("," , map {"snpeff_$_"} @header_col) . " " ;
#$cmd_gannot .= " -t " . join("," , map {"text"} @header_col) . " " ;
#$cmd_gannot .= " -c " . join("," , map {"snpeff_$_"} @header_col) . " " ;
#$cmd_gannot .= " -o " . join("," , map {"first"} @header_col) . " " ;
#$cmd_gannot .= " $file_geminidb 2> $dir_log/stderr_gannot_snpeff_$jobid.log" ;
#$log->andRun($cmd_gannot) ;
#

# Step 5. alt. store list of annotation column and run gemini-annotate later
my $output_list = {
    extract_list   => [map {"snpeff_$_"} @header_col ] ,
    type_list      => [map {"text"}      @header_col ] ,
    coladd_list    => [map {"snpeff_$_"} @header_col ] ,
    operation_list => [map {"first"}      @header_col] ,
} ;

store $output_list , $file_hashref_vannot_tmp ; 

