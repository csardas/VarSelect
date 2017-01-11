#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use Storable ;
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use VarSelect::Vcf_parser ;

my $DEFAULT_CHUNKSIZE = 10000 ;
#my $DEFAULT_CHUNKSIZE = 100000 ;

# default global setting
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

# options handle
my $opts = {} ;
getopts("j:d:v:g:h", $opts) ;

my $jobid		     = $opts->{j} ;
my $file_db		     = $opts->{d} ;
my $file_vcfgz_input	     = $opts->{v} ;
my $file_hashref_gannot_cols = $opts->{g} ;

my $dir_results_output  = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log		= "$dir_results_output/log_$jobid" ;
my $dir_work		= "$dir_results_output/work_$jobid" ;

my $file_log = "$dir_log/log_gannot_vsannot_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("Gemini annotate from annotation of vs_annot start") ;


my $gannot_cols = retrieve $file_hashref_gannot_cols ;
# extract_list
# type_list
# coladd_list
# operation_list


# Step 0. get var_num
my $var_num = `zgrep -c -v '^#' $file_vcfgz_input ` ;
chomp $var_num ;
print "var num: $var_num\n" ;
print "Default chunk size: $DEFAULT_CHUNKSIZE \n" ;

if ($var_num <=  $DEFAULT_CHUNKSIZE) {
    my $cmd_gannot = "gemini annotate " ;
    $cmd_gannot .= " -f $file_vcfgz_input " ;
    $cmd_gannot .= " -a extract " ;
    $cmd_gannot .= " -e " . join("," , @{$gannot_cols->{extract_list}}) ;
    $cmd_gannot .= " -t " . join("," , @{$gannot_cols->{type_list}}) ;
    $cmd_gannot .= " -c " . join("," , @{$gannot_cols->{coladd_list}}) ;
    $cmd_gannot .= " -o " . join("," , @{$gannot_cols->{operation_list}}) ;
    $cmd_gannot .= " $file_db " ;
    $cmd_gannot .= " 2> $dir_log/stderr_gannot_vsanno_$jobid.log" ;

    $log->andRun($cmd_gannot) ;

} else {  # if varnum > DEFAULT_CHUNKSIZE
    # Step 1. extract header
    my $file_header = "$dir_work/vcf_header_for_gannot1.vcf" ;
    my $cmd_create_header = "zgrep '^#' $file_vcfgz_input > $file_header"  ;
    $log->andRun($cmd_create_header) ;

    my $chunk_id = 0 ;
    my $input_number = 0 ;
    my $input_number_all = 0 ;
    
    my $file_input_chunk = "$dir_work/gannot_vcf_input_$chunk_id.vcf" ;

    open (my $TGT_vcf_input_chunk,">$file_input_chunk" ) ;

    open (my $SRC_vcfinput,"zgrep -v '^#' $file_vcfgz_input|") ;
    while (my $line = <$SRC_vcfinput>) {
	$input_number ++ ;
	$input_number_all ++ ;
	print $TGT_vcf_input_chunk $line ;

	if ($input_number == $DEFAULT_CHUNKSIZE || eof($SRC_vcfinput) ) {
	    close $TGT_vcf_input_chunk ;
	    $input_number = 0 ;

	    my $cmd_combine_header_vcf = "cat $file_header $file_input_chunk | bgzip -@ 32 -c > $file_input_chunk.gz " ;
	    $log->andRun($cmd_combine_header_vcf) ;
	    tabix_vcf("$file_input_chunk.gz") ;

	    my $cmd_gannot = "gemini annotate " ;
	    $cmd_gannot .= " -f $file_input_chunk.gz " ;
	    $cmd_gannot .= " -a extract " ;
	    $cmd_gannot .= " -e " . join("," , @{$gannot_cols->{extract_list}}) ;
	    $cmd_gannot .= " -t " . join("," , @{$gannot_cols->{type_list}}) ;
	    $cmd_gannot .= " -c " . join("," , @{$gannot_cols->{coladd_list}}) ;
	    $cmd_gannot .= " -o " . join("," , @{$gannot_cols->{operation_list}}) ;
	    $cmd_gannot .= " $file_db " ;
	    $cmd_gannot .= " 2> $dir_log/stderr_gannot_vsanno_$jobid.chunk$chunk_id.log" ;
	    $log->andRun($cmd_gannot) ;

	    print "update $input_number_all vars\n" ;

	    $chunk_id ++ ;
	    $file_input_chunk = "$dir_work/gannot_vcf_input_$chunk_id.vcf" ;
	    open ($TGT_vcf_input_chunk , ">$file_input_chunk") ;
	}
    }
    
}

#
#
# loop
# extract first chunk
# ensembl chunk_x.vcf
# gannot

$log->write("Gemini annotate from annotation of vs_annot finish.") ;
