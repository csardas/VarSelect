#!/usr/bin/perl
use strict ;
use Getopt::Std ;
use File::Basename ;
use File::Path qw(make_path remove_tree);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0) . '/lib';
use VarSelect ;
use VarSelect::Log ;

# should be put in global config
my $DEFAULT_threads = 16 ;

my $jobid = getts() ;

# Command handle
my $commands_avail = {map {$_ => 1} qw/annotate analysis compare/}  ;
my $command = shift @ARGV ;
die "\nOnly " . join (", ", keys %$commands_avail ) . " commands are allowed!\n\n" . usage_all() unless exists $commands_avail->{$command} ;

my $workflow_avail = {map {$_ => 1} qw/family paired none/} ;

# VarSelect Path settings
my $dir_vsHOME = dirname(abs_path($0)) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_scripts = "$dir_vsHOME/scripts" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

my $script_vsannotate = "$dir_scripts/vs_annotate.pl" ;
my $script_vsanalysis = "$dir_scripts/vs_analysis.pl" ;
my $script_vscompare  = "$dir_scripts/vs_compare.pl" ;

my $opts = {} ;
if ($command eq 'annotate') {

    # Options for annotate
    getopts("v:p:m:c:x:hkuin:" , $opts) ;

    die usage_annotate()  if ($opts->{h}) ;
    my $file_vcflist	    = $opts->{v} || die "-v vcf_file_list is required!\n\n" . usage_annotate() ;
    my $file_ped	    = $opts->{p} || die "-p ped file is required!\n\n" . usage_annotate() ;
    my $mode_workflows	    = $opts->{m} || die "-m analysis mode is required, [paired, family, none] !\n\n" . usage_annotate()  ;
    my $flag_multicallers   = $opts->{k} ;
    my $threads_enable	    = $opts->{n} || $DEFAULT_threads ;
    my $file_cnv	    = $opts->{c} ;
    my $file_xpr	    = $opts->{x} ;
#    my $freq_vgt_filter	    = $opts->{f} || 1 ;

    my $flag_combine_union  = $opts->{u} ;
    my $flag_combine_intersect = $opts->{i} ;

    my $type_combine = "union" ;
    if ($flag_combine_union && $flag_combine_intersect) {
	die "\nError: -u and -i can NOT enable at same time!\n\n" . usage() ;
    } elsif ($flag_combine_intersect) {
	$type_combine = "intersection" ;
    }

    die "Only [ " . join (",", sort keys %$workflow_avail) . " ] are available for option -m \n" unless (exists $workflow_avail->{$mode_workflows}) ;

    # Output setting
    my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
    my $dir_log = "$dir_results_output/log_$jobid" ;
    my $dir_work = "$dir_results_output/work_$jobid" ;

    die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
    make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
    make_path($dir_log , $dir_work ,{ chmod => 0755 ,} )  ;

    my $file_stderr = "$dir_log/stderr_vs_annotate.$jobid.log" ;
    my $file_log = "$dir_log/log_varselect_$jobid.log" ;
    my $log = VarSelect::Log->new(file => $file_log) ;
    $log->write("VarSelect Start") ;
    $log->write("Jobid: $jobid") ;
    $log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ) ;

    my ($fname_vfile,$fpath_vfile,$fext_vfile) = fileparse($file_vcflist , qw/.txt/) ;
    my $file_db = "./$fname_vfile\_varselect.db" ;


    print "VarSelect Job id: $jobid\n" ;

    # Step 1. Annotate from vcf
    my $cmd_vs_annotate = "$script_vsannotate " ;
    $cmd_vs_annotate .= " -v $file_vcflist " ;
    $cmd_vs_annotate .= " -p $file_ped " ;
    $cmd_vs_annotate .= " -j $jobid " ;
    $cmd_vs_annotate .= " -n $threads_enable " ;
    $cmd_vs_annotate .= " -k " if ($flag_multicallers) ;
    $cmd_vs_annotate .= " -d $file_db " ; # make sure where to put db file
    $cmd_vs_annotate .= " 2> $file_stderr " ;

    $log->write("VS_annotate start") ;
    $log->andRun($cmd_vs_annotate) ;
    $log->loadlog("$dir_log/log_vsannotate_$jobid.log") ;
    $log->write("VS_annotate finish") ;

    # Step 2. Analysis according to PED info
    my $cmd_vs_analysis = "$script_vsanalysis " ; 
    $cmd_vs_analysis .= " -j $jobid " ;
    $cmd_vs_analysis .= " -d $file_db " ;
    $cmd_vs_analysis .= " -p $file_ped " ;
    $cmd_vs_analysis .= " -m $mode_workflows " ;
#    $cmd_vs_analysis .= " -f $freq_vgt_filter " ;
    $cmd_vs_analysis .= " -k " if ($flag_multicallers) ;
    $cmd_vs_analysis .= " -i " if ($flag_combine_intersect) ;
    $cmd_vs_analysis .= " -u " if ($flag_combine_union) ;
    $cmd_vs_analysis .= " -c $file_cnv " if ($file_cnv) ;
    $cmd_vs_analysis .= " -x $file_xpr " if ($file_xpr) ;

    $log->write("VS_analysis start") ;
    $log->andRun($cmd_vs_analysis) ;
#    $log->loadlog("$dir_log/log_vsanalysis_$jobid.log") ;
    $log->write("VS_analysis finish") ;

} elsif ($command eq 'analysis') {
    # Options for analysis
    getopts("d:p:m:hkuic:x:n:f:",$opts) ;

    die usage_analysis() if $opts->{h} ;

    my $file_db		= $opts->{d} || die "\n\n-d is required!\n" . usage_analysis() ;
    my $file_ped	= $opts->{p} || die "\n\n-p is required!\n" . usage_analysis() ;
    my $mode_workflows	= $opts->{m} || die "\n\n-m is required!\n" . usage_analysis() ;

    my $flag_multicallers   = $opts->{k} ;
    my $threads_enable	    = $opts->{n} || $DEFAULT_threads ;
    my $file_cnv	    = $opts->{c} ;
    my $file_xpr	    = $opts->{x} ;
    my $freq_vgt_filter	    = $opts->{f} || 1;

    my $flag_combine_union  = $opts->{u} ;
    my $flag_combine_intersect = $opts->{i} ;

    my $type_combine = "union" ;
    if ($flag_combine_union && $flag_combine_intersect) {
	die "\nError: -u and -i can NOT enable at same time!\n\n" . usage() ;
    } elsif ($flag_combine_intersect) {
	$type_combine = "intersection" ;
    }

    die "\ndb $file_db doesn't exist!\n" unless (-e $file_db) ;
    die "\ndb $file_db is not writable!\n" unless (-w $file_db) ;

    die "Only [ " . join ("," , sort keys %$workflow_avail) . " ] are available for option -m \n" unless (exists $workflow_avail->{$mode_workflows}) ;

    print "VarSelect Job id: $jobid\n" ;

    # Output setting, create output directory
    my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
    my $dir_log = "$dir_results_output/log_$jobid" ;
    my $dir_work = "$dir_results_output/work_$jobid" ;

    die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
    make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
    make_path($dir_log , $dir_work ,{ chmod => 0755 ,} )  ;

    print "$dir_results_output is created.\n" ;

    my $file_stderr = "$dir_log/stderr.vs_analysis.$jobid.log" ;

    # Start logging
    my $file_log = "$dir_log/log_varselect_$jobid.log" ;
    my $log = VarSelect::Log->new(file => $file_log) ;
    $log->write("Jobid: $jobid") ;
    $log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ) ;

    # pipeline setting
    my $cmd_vs_analysis = "$script_vsanalysis " ; 
    $cmd_vs_analysis .= " -j $jobid " ;
    $cmd_vs_analysis .= " -d $file_db " ;
    $cmd_vs_analysis .= " -p $file_ped " ;
    $cmd_vs_analysis .= " -m $mode_workflows " ;
#    $cmd_vs_analysis .= " -f $freq_vgt_filter " ;
    $cmd_vs_analysis .= " -k " if ($flag_multicallers) ;
    $cmd_vs_analysis .= " -i " if ($flag_combine_intersect) ;
    $cmd_vs_analysis .= " -u " if ($flag_combine_union) ;
    $cmd_vs_analysis .= " -c $file_cnv " if ($file_cnv) ;
    $cmd_vs_analysis .= " -x $file_xpr " if ($file_xpr) ;

    $log->write("VS_analysis start") ;
    $log->andRun($cmd_vs_analysis) ;
    $log->loadlog("$dir_log/log_vsanalysis_$jobid.log") ;
    $log->write("VS_analysis finish") ;

    print "VarSelect $jobid done.\n" ;
    print localtime() ."\n" ;

} elsif ($command eq 'compare') {
    getopts("a:b:c:h",$opts) ;

    die usage_compare() if ($opts->{h}) ;
    my $datasetA	= $opts->{a} || die "\n-a is required!\n" . usage_compare() ;
    my $datasetB	= $opts->{b} || die "\n-b is required!\n" . usage_compare() ;
    my $method_compare	= $opts->{c} || die "\n-c is required!\n" . usage_compare() ;

    # Output setting
    my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
    my $dir_log = "$dir_results_output/log_$jobid" ;
    my $dir_work = "$dir_results_output/work_$jobid" ;

    die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
    make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
    make_path($dir_log , $dir_work ,{ chmod => 0755 ,} )  ;

    my $file_stderr = "$dir_log/stderr_vs_compare.$jobid.log" ;
    my $file_log = "$dir_log/log_varselect_$jobid.log" ;
    my $log = VarSelect::Log->new(file => $file_log) ;
    $log->write("VarSelect Start") ;
    $log->write("Jobid: $jobid") ;
    $log->write("VS compare start") ;
    $log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ) ;

    my $cmd_vs_compare = "$script_vscompare " ;
    $cmd_vs_compare .= " -a $datasetA " ;
    $cmd_vs_compare .= " -b $datasetB " ;
    $cmd_vs_compare .= " -c $method_compare " ;
    $cmd_vs_compare .= " -j $jobid " ;
    
    $log->andRun($cmd_vs_compare) ;
    $log->write("VS compare finish") ;

} elsif ($command eq 'joblist') {
    getopts("d:",$opts) ;
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

    # PED
    # ID
    # columns
    return $output ;
}
    

sub usage_all {
return <<EOusage
$0 <command> [options]

Commands:
    annotate - Annotate VCFs from scratch
    analysis - Re-analysis with different settings
    compare  - Compare results between analysis

use $0 <command> -h to get more detail help page
EOusage
}

sub usage_annotate {
    return <<EOusage1
$0 annotate 
    -v VCF files list					 (required)
    -p PED file						 (required)
    -m workflow mode: [paired, family, none]		 (required)

    -k multiple caller mode				 (optional)
    -u get union sets of variants between callers	 (only work with -k option)
    -i get intersection sets of variants between callers (only work with -k option)

    -c CNVkit log2 calls cns file list			 (optional)
    -x expression profiles list				 (optional)
EOusage1
}

sub usage_analysis {
    return <<EOusage2
$0 analysis

    -d db file						 (required)
    -p PED file						 (required)
    -m workflow mode: [paired, family, none]	    	 (required)

    -k multiple caller mode				 (optional)
    -u get union sets of variants between callers	 (only work with -k option)
    -i get intersection sets of variants between callers (only work with -k option)

    -c CNVkit log2 calls cns file list			 (optional)
    -x expression profiles list				 (optional)
    
EOusage2
}

sub usage_compare {
    return <<EOusage3
$0 compare 

    -a Path/to/result/of/datasetA
    -b Path/to/result/of/datasetB

    -c compare mode [1-4]
       1: A or B
       2: A and B
       3: A only
       4: B only
EOusage3
}

