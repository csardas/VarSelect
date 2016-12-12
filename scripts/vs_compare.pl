#!/usr/bin/perl
use strict ;

use Getopt::Std ;
use File::Basename ;
use Storable ;
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;

my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

my $method_compare_avail= {
    1 => 'A and B' ,
    2 => 'A or B' ,
    3 => 'A only' ,
    4 => 'B only' ,
} ;

my $opts = {} ;

    getopts("a:b:c:j:",$opts) ;

    my $jobid = $opts->{j} ;
    my $dir_setA        = $opts->{a} || die "\n-a is required\n\n". usage_compare() ;
    my $dir_setB        = $opts->{b} || die "\n-b is required\n\n". usage_compare() ;
    my $method_compare    = $opts->{c} ;

    die "$dir_setA is not exists!\n" . usage_compare() unless (-e $dir_setA) ;
    die "$dir_setB is not exists!\n" . usage_compare() unless (-e $dir_setB) ;
    die "-c should be 1-4!\n" . usage_compare() unless (exists $method_compare_avail->{$method_compare}) ;

    my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
    my $dir_log = "$dir_results_output/log_$jobid" ;
    my $dir_work = "$dir_results_output/work_$jobid" ;

    my $file_log = "$dir_log/log_vscompare_$jobid.log" ;
    my $log = VarSelect::Log->new(file => $file_log) ;

    my $DEFAULT_varselect_config = "varselect.config" ;

    my $file_config_a = "$dir_setA/$DEFAULT_varselect_config" ;
    my $file_config_b = "$dir_setB/$DEFAULT_varselect_config" ;
    my $configA = load_config($file_config_a) ;
    my $configB = load_config($file_config_b) ;

    die "dataset A and B come from different DB!\n\tDB for datasetA: $configA->{DB}\n\tDB for datasetB: $configB->{DB}\n" if ($configA->{DB} ne $configB->{DB}) ;
    my $gemini_db = $configA->{DB} ;
    print "gemini_db = $gemini_db\n" ;
    my ($dbfn,$dbpath,$dbext) = fileparse($gemini_db , qw/.db/) ;
    my $file_vcfgz_totalannot = $dbpath . $dbfn . ".vcf.gz" ;
    my $file_vcf_finaloutput = "VarSelect_compare_$jobid.vcf" ;
    my $file_txt_finaloutput = "VarSelect_compare_$jobid.txt" ;

    die "$gemini_db should have coresspond vcf $file_vcfgz_totalannot ." unless (-e $file_vcfgz_totalannot) ;
    my $vsidA = $configA->{ID} ;
    my $vsidB = $configB->{ID} ;

    my $iaA = "in_analysis_$vsidA" ;
    my $iaB = "in_analysis_$vsidB" ;

    my $compare_mode = '' ;
    my $query_sql = "select chrom,start+1,ref,alt,1 from variants where " ;
    my $query_sql_dump = "select * from variants where " ;
    if ($method_compare == 1) {      
        $query_sql .= " $iaA = \"1\" or $iaB = \"1\"" ;
        $query_sql_dump .= " $iaA = \"1\" or $iaB = \"1\"" ;
	$compare_mode = 'union' ;

    } elsif ($method_compare == 2) {
        $query_sql .= " $iaA = \"1\" and $iaB = \"1\"" ;
        $query_sql_dump .= " $iaA = \"1\" and $iaB = \"1\"" ;
	$compare_mode = 'intersection' ;

    } elsif ($method_compare == 3) {
        $query_sql .= " $iaA = \"1\" and $iaB = \"0\"" ;
        $query_sql_dump .= " $iaA = \"1\" and $iaB = \"0\"" ;
	$compare_mode = 'Aonly' ;

    } elsif ($method_compare == 4) {
        $query_sql .= " $iaA = \"0\" and $iaB = \"1\"" ;
        $query_sql_dump .= " $iaA = \"0\" and $iaB = \"1\"" ;
	$compare_mode = 'Bonly' ;

    } else {  # should not be here
    }
    
    my $result_compare_tab_annot = "$dir_results_output/vs_compare_output_annot_$jobid.txt" ;
    my $result_compare_tab = "$dir_results_output/vs_compare_output_$jobid.txt" ;
    my $result_compare_vcf = "$dir_results_output/vs_compare_output_$jobid.vcf" ;

    gemini_query ("$query_sql" , $gemini_db , "$result_compare_tab_annot" ,0) ;
    gemini_query ("$query_sql_dump" , $gemini_db , "$result_compare_tab" ,1) ;
    gemini_query_vcf ("$query_sql_dump" , $gemini_db , $result_compare_vcf , 1) ;

#    my $cmd_bgzip = "bgzip $result_compare_vcf -c > $result_compare_vcf.gz"  ;
#    tabix_vcf($result_compare_tab ,$log) ;
    tabix_vcf($result_compare_tab_annot ,$log) ;
    tabix_vcf($result_compare_vcf ,$log) ;

    my $file_vcfgz_outputtmp = "$dir_work/compare_output_tmp.vcf.gz" ;
    my $cmd_vannot = "cat $result_compare_vcf | vcf-annotate " ;
    $cmd_vannot .= " -a $result_compare_tab_annot.gz " ;
    $cmd_vannot .= " -c " . join ("," , (qw/CHROM POS REF ALT/ ,"INFO/in_analysis_$jobid") ) ;
    $cmd_vannot .= " -d  key=INFO,ID=in_analysis_$jobid,Number=1,Type=Integer,Description='Secondary analysis $jobid ' " ;
    $cmd_vannot .= " 2> $dir_log/stderr_vannot_compare.$jobid.log | bgzip -c > $file_vcfgz_outputtmp" ;
    $log->andRun($cmd_vannot) ;

    tabix_vcf ($file_vcfgz_outputtmp ,$log) ;

    my $cmd_gannot = "gemini annotate " ;
    $cmd_gannot .= " -f $file_vcfgz_outputtmp " ;
    $cmd_gannot .= " -a extract " ;
    $cmd_gannot .= " -c in_analysis_$jobid "  ;
    $cmd_gannot .= " -e in_analysis_$jobid " ;
    $cmd_gannot .= " -t integer " ;
    $cmd_gannot .= " -o first " ;
    $cmd_gannot .= " $gemini_db " ;
    $cmd_gannot .= " 2> $dir_log/stderr_gannot_$jobid.log" ;
    $log->andRun($cmd_gannot) ;

sub gemini_query_vcf {
    my $sql = shift ;
    my $db = shift ;
    my $file_output = shift ;
    my $header = shift || 0 ;

    my $cmd_gemini_query = "gemini query " ;
    $cmd_gemini_query .= " --format vcf " ;
    $cmd_gemini_query .= " --header " if ($header) ;
    $cmd_gemini_query .= " -q '$sql' $db > $file_output" ;
    $log->andRun($cmd_gemini_query) ;
}


sub gemini_query {
    my $sql = shift ;
    my $db = shift ;
    my $file_output = shift ;
    my $header = shift || 0 ;

    my $cmd_gemini_query = "gemini query " ;
    $cmd_gemini_query .= " --header " if ($header) ;
    $cmd_gemini_query .= " -q '$sql' $db > $file_output" ;
    $log->andRun($cmd_gemini_query) ;
}


sub load_config {
    my $file = shift ;
    die "$file isn't exist! please check it's a varselect_analysis_dir or not" unless (-e $file) ;

    my $output = {} ;
    open (my $SRC , "$file") ;
    while (<$SRC>) {
        chomp ;
        my ($key,$val) = split (/\t/ , $_ , 2) ;
        if (exists $output->{$key}) {
            $output->{$key} .= "\n$val" ;
        } else {
            $output->{$key} = $val ;
        }
    }
    close $SRC ;

    return $output ;
}
    
