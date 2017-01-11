#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use File::Path qw(make_path remove_tree);
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

my $ts = getts() ;
my $DEFAULT_build = $Setting->{genome_annovar} ;
my $dir_annovar = $Setting->{dir_annovar} ;
my $exe_annovar = "$dir_annovar/table_annovar.pl" ;
my $file_humandb_annovar = "$dir_annovar/humandb" ;


# Options handling
my $opts = {} ;
getopts("i:o:l:j:", $opts) ;

my $jobid = $opts->{j} ;
my $file_input = $opts->{i} ;
my $file_output = $opts->{o} ;

my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;
my $file_config = "$dir_results_output/varselect.config" ;
my $file_log = "$dir_log/running_annovar_$jobid.log" ;

my $file_prefix = "$dir_results_output/runannovar_$ts" ;

my $log = VarSelect::Log->new(file => $file_log) ;

$log->write ("Running ANNOVAR start") ;
my $file_input_annovar = $file_input ;
my $file_vcf_annovar = "$file_prefix.$DEFAULT_build\_multianno.vcf" ;

my $protocal_enable = {
    avsift => 'f' ,
    dbscsnv11 => 'f' ,
    avsnp147 => 'f' ,
    cg69 => 'f' ,
    cosmic70 => 'f' ,
    esp6500siv2_all => 'f' ,
    exac03 => 'f' ,
    exac03nontcga => 'f' ,
    exac03nonpsych => 'f' ,
    kaviar_20150923 => 'f' ,
    hrcr1 => 'f' ,
    '1000g2015aug_all' => 'f' ,
    gme => 'f' ,
    nci60 => 'f' ,
    icgc21 => 'f' ,
    clinvar_20160302 => 'f' ,
    popfreq_max_20150413 => 'f' ,
    popfreq_all_20150413 => 'f' ,
    mitimpact24 => 'f' ,
    'gerp++gt2' => 'f' ,
    'phastConsElements46way' => 'r' ,
    'dbscsnv11' => 'f' ,
    'tfbsConsSites' => 'r' ,
    'gwasCatalog' => 'r' ,
    'wgRna' => 'r' ,
    fathmm => 'f' ,

} ;

my $list_protocols = [] ;
my $list_operation = [] ;

foreach my $protocal (keys %$protocal_enable) {
    push @$list_protocols , $protocal ;
    push @$list_operation , $protocal_enable->{$protocal} ;
}

my $cmd_annovar = "$exe_annovar $file_input_annovar $file_humandb_annovar " ;
$cmd_annovar .= " -buildver $DEFAULT_build " ;
$cmd_annovar .= " -out $file_prefix " ;
$cmd_annovar .= " -remove " ;
$cmd_annovar .= " -protocol " . join ("," , map {"'" . $_ . "'" } @$list_protocols ) ;
$cmd_annovar .= " -operation " . join("," , @$list_operation ) ;
# $cmd_annovar .= " -nastring . " ;
$cmd_annovar .= " -vcfinput " ;
$cmd_annovar .= " 2> $dir_log/stderr.annovar.$ts.log" ;

$log->write ("$cmd_annovar") ;
`$cmd_annovar` ;
$log->write ("Running ANNOVAR finished") ;

my $cmd_bgzip_output = "bgzip -c $file_prefix.$DEFAULT_build\_multianno.vcf > $file_output" ;
$log->write("$cmd_bgzip_output") ;
`$cmd_bgzip_output` ;




