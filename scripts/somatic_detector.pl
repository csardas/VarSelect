#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use VarSelect::Vcf_parser ;

# default global setting
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

# Options handle
my $opts = {} ;
getopts("i:o:p:j:d:n:", $opts) ;

my $jobid		  = $opts->{j} || getts() ;
my $file_db		  = $opts->{d} ;
my $file_ped		  = $opts->{p} || die "\nError: -p PED file is required!\n\n" ;
my $file_vcf_input	  = $opts->{i} || die "\nError: -i VCF file is required!\n\n" ;
my $file_output_isSomatic = $opts->{o} ;
my $threads_enable	  = $opts->{n} ;

my ($input_name,$input_path,$input_ext) = fileparse($file_vcf_input, qw/.vcf .vcf.gz/) ;
my ($fname_db,$fpath_db,$fext_db) = fileparse($file_db, qw/.db/) ;

# jobs specific path
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
make_path($dir_log , $dir_work ,{ chmod => 0755 ,} ) unless (-e $dir_work && -e $dir_log) ;

my $file_vcfgz_isSomatic = "$dir_results_output/$input_name\_somatic\_$jobid.vcf.gz" ;

# start logging
my $file_log = "$dir_log/log_somaticdetect_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("Somatic detecting start") ;

# Step 1. Load PED sample info
#=====================================
open (my $SRC_ped, "$file_ped") || die "\nCan't open PED file $file_ped !\n\n" ;

my $list_affected_samples = [] ;
my $list_unaffected_samples = [] ;
while (<$SRC_ped>) {
    chomp ;
    next if /^#/ ;
    my ($family,$sample,$father,$mother,$sex,$phenotype) = split /\t/ ;
    push @$list_affected_samples , $sample if ($phenotype == 2) ;
    push @$list_unaffected_samples , $sample if ($phenotype == 1) ;
}

close $SRC_ped ;

# Step 2. Load VCF 
#=====================================
my $vcf_parser = Vcf_parser->new (file => $file_vcf_input) ;
open (my $TGT_rst , ">$file_output_isSomatic") ;
print $TGT_rst join("\t", (qw/#CHROM POS REF ALT is_Somatic/ , @{$vcf_parser->{samples}}) ) . "\n"  ;

while (my $var = $vcf_parser->next_var) {
    # Step 3. list all un-affected genotype
    my $is_Somatic = 1 ;
    my $record_uaf_gt = {} ;
    foreach my $sample_uaf (@$list_unaffected_samples) {
	my $gt_uaf = gtchk0($var->{sample_val}->{$sample_uaf}->{GT}) ;
	$record_uaf_gt->{$gt_uaf} ++ ;
    }

    # Step 4. check each affected sample, drop while equal any u-af gt
    my $af_gt_list = {} ;
    foreach my $sample_af (@$list_affected_samples) {
	my $gt_af = gtchk0($var->{sample_val}->{$sample_af}->{GT}) ;

	if ($gt_af eq '0/0' || exists $record_uaf_gt->{$gt_af}) {
	    $is_Somatic = 0 ;
	} else {
	    # somatic
	    $af_gt_list->{$gt_af} ++ ;
	}
    }

    my $somatic_gt = '-' ;

    if ($is_Somatic) {
	$somatic_gt = join("," , sort keys %$af_gt_list ) ;
    }

    print $TGT_rst join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) ;
    print $TGT_rst "\t$is_Somatic\t$somatic_gt" ;

    foreach my $sample (@{$var->{samples}}) {
        print $TGT_rst "\t$sample:" . $var->{sample_val}->{$sample}->{GT} ;
    }

    print $TGT_rst "\n" ;
}

tabix_vcf($file_output_isSomatic) ;

# Step 5. vcf-annotate 
my $cmd_vannot = "zcat $file_vcf_input | vcf-annotate -a $file_output_isSomatic.gz " ;
$cmd_vannot .= " -c " . join("," , qw"CHROM POS REF ALT INFO/is_Somatic INFO/Somatic_gt" ) ;
$cmd_vannot .= " -d key=INFO,ID=is_Somatic,Number=1,Type=Integer,Description='is_Somatic annotation by Varselect' " ;
$cmd_vannot .= " -d key=INFO,ID=Somatic_gt,Number=.,Type=String,Description='Somatic_gt annotation by Varselect' " ;
$cmd_vannot .= " 2> $dir_log/stderr_vannot_somatic_$jobid.log " ;
$cmd_vannot .= " | bgzip -@ $threads_enable -c > $file_vcfgz_isSomatic " ;

$log->andRun($cmd_vannot) ;
tabix_vcf($file_vcfgz_isSomatic) ;

# Step 6. gemini annotate
my $cmd_gannot = "gemini annotate -f $file_vcfgz_isSomatic " ;
$cmd_gannot .= " -a extract " ;
$cmd_gannot .= " -e is_Somatic,Somatic_gt " ;
$cmd_gannot .= " -t integer,text " ;
$cmd_gannot .= " -c is_Somatic_$jobid,Somatic_gt_$jobid" ;
$cmd_gannot .= " -o first,first " ;
$cmd_gannot .= " $file_db 2> $dir_log/stderr_gannot_somatic_$jobid.log" ;
$log->andRun($cmd_gannot) ;

$log->write("Somatic detecting finish") ;


sub gtchk0 {
    my $gt = shift ;
    $gt = '0/0' if ($gt eq '.' || !defined($gt) ) ;
    $gt =~ s/\./0/g ; # treat . as 0
    return $gt ;
}

