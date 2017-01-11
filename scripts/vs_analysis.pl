#!/usr/bin/perl
use strict ;
use Getopt::Std ;
use Storable ;
use File::Basename ;
use File::Path qw(make_path) ;
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

# options handle
my $opts = {} ;
getopts("j:d:p:m:n:kuix:c:",$opts) ;
my $jobid               = $opts->{j} ;
my $file_db		= $opts->{d} ;
my $file_ped            = $opts->{p} ;
my $mode_workflow	= $opts->{m} ;
my $flag_multicallers   = $opts->{k} || 0 ;
my $threads_enable      = $opts->{n} || 16;
my $file_cnv_list	= $opts->{c} ;
my $file_xpr_list	= $opts->{x} ;

my $type_combine = "union" ;
if ($opts->{u} && $opts->{i}) {
    die "\nError: -u and -i can NOT enable at same time!\n\n" . usage() ;
} elsif ($opts->{i}) {
    $type_combine = "intersection" ;
}

# jobs specific path
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
make_path($dir_log , $dir_work ,{ chmod => 0755 ,} ) unless (-e $dir_work && -e $dir_log) ;

my $file_config = "$dir_results_output/varselect.config" ;

my ($fname_db, $fpath_db , $file_ext) = fileparse($file_db, qw/.db/) ;
my $file_prefix = "$dir_results_output/$fname_db" ;
my $file_config_db = "$fpath_db/$fname_db.config" ;

# Start logging
my $file_log = "$dir_log/log_vsanalysis_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("VSanlz start") ;
$log->write("Jobid: $jobid") ;
$log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ) ; # list of options

my $config_from_db = {} ;
$log->write("Loading db config file: $file_config_db") ;
open (my $SRC_config , "$file_config_db" ) || die "cant loading $file_config_db" ;
while (<$SRC_config>) {
    chomp ;
    my($idx,$val) = split /\t/ ;
    $config_from_db->{$idx} = $val ;
}
close $SRC_config ;

my $file_hashref_gannot_tmp = "$dir_work/gannot_tmp_vsanal.hashref" ;

open (my $CONFIG , ">$file_config") ;
print $CONFIG "ID\t$jobid\n" ;
print $CONFIG "DB\t$file_db\n" ;
print $CONFIG "caller_list\t$config_from_db->{caller_list}\n" if (exists $config_from_db->{caller_list}) ;
print $CONFIG "sample_list\t$config_from_db->{sample_list}\n" if (exists $config_from_db->{sample_list}) ;

my $caller_list = {map {$_ => 1} split(/\,/,$config_from_db->{caller_list}) } ;
my $sample_list = {map {$_ => 1} split(/\,/,$config_from_db->{sample_list}) } ;

# Update PED
my $file_sample_sex = $file_prefix ."_sample_sex.txt" ;
my $file_sample_sex_mcaller = $file_prefix ."_sample_sex_mcaller.txt" ;
my $list_samples = {} ;

my $cmd_gemini_amend = "gemini amend --sample $file_ped $file_db 2> $dir_log/stderr_ped_amend_$jobid.log" ;
$log->andRun ($cmd_gemini_amend) ;

my ($ped_fn,$ped_path,$ped_ext) = fileparse($file_ped,qw/.ped/) ;
my $cmd_cp_ped_to_anadir = "cp $file_ped $dir_results_output/$ped_fn\_$jobid.ped" ;

# Loading PED
# check affected status of sample
my $affect_samples = [] ;
my $unaffect_samples = [] ;
my $samples = [] ;
open (my $SRC_ped,"$file_ped") ;
while (my $line = <$SRC_ped>) {
    chomp $line;
    my ($family,$sample,$father,$mother,$sex,$affect_status) = split(/\t/,$line) ;

    if ($affect_status == 1) {
	push @$unaffect_samples , $sample ;

    } elsif ($affect_status == 2) {
	push @$affect_samples , $sample ;
    }

    push @$samples , $sample ;
}
close $SRC_ped ;

# Dump VCF from gemini db 
my $file_vcfgz_foranalysis = "$fpath_db/$fname_db.vcf.gz" ; # default setting

if (-e $file_vcfgz_foranalysis) {
    tabix_vcf($file_vcfgz_foranalysis,$log) ; # unless (-e "$file_vcfgz_foranalysis.tbi") ;

} else {
    die("$file_vcfgz_foranalysis not found! please check your db create directory and put the $fname_db.vcf.gz to the current location of your db!\n") ;
}



# Check multi-caller sample
if ($flag_multicallers) {
    if ($type_combine eq "union") {

        my $vcf_parser = Vcf_parser->new(file=> $file_vcfgz_foranalysis) ;

        my $file_output_union_conflict = "$dir_results_output/multicaller_union_inconsistant_$jobid.txt" ;
        my $file_vcf_union_set = "$dir_results_output/multicaller_unionset_$jobid.vcf";

        open (my $TGT_union_conflict , ">$file_output_union_conflict" ) ;
        open (my $TGT_union_vcf, ">$file_vcf_union_set") ;

        print $TGT_union_vcf join("\t" , (qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/ , sort keys %$sample_list)) . "\n" ;

        while (my $var = $vcf_parser->next_var) {
            my $flag_inanalysis = 1 ;
            my $reason_drop = '' ;
            my $gt_all_in_one_locus = [] ;
            my $list_gt_by_sample = [] ;

           # conflict if more than one gt type in gt_got
           foreach my $sample (sort keys %$sample_list) {
                my $gt_got = {} ;
                foreach my $caller (sort keys %$caller_list) {
                    my $sample_tag = $sample . "_" . $caller ;
                    my $gt = $var->{sample_val}->{$sample_tag}->{GT} ;
                    $gt =~ s/\|/\// ;   # ignore phased, treat as an un-phased genotype
                    next if ($gt eq '0/0') ;
                    $gt_got->{$gt} ++ ;

                    push @$gt_all_in_one_locus , "$sample_tag:$gt" ;
                }

                my $list_of_got_gt = [keys %$gt_got] ;
                my $gt_type_count = scalar @$list_of_got_gt ;

                if ($gt_type_count > 1) {
		    #union con-flict
                   $flag_inanalysis = 0 ;
                    $reason_drop .= "Conflict genotype in 1+ callers of sample $sample,"

                } elsif ($gt_type_count == 1) {
                    push @$list_gt_by_sample , $list_of_got_gt->[0] ;

                } else {
                    push @$list_gt_by_sample , '0/0' ;
                }
            }

            if ($flag_inanalysis) {
                # output vcf with sample only (no caller now) for workflow analysis
                print $TGT_union_vcf join("\t" , map {$var->{$_} } qw/CHROM POS ID REF ALT QUAL FILTER/ ) . "\t" . "\tGT\t" . join("\t",@$list_gt_by_sample) . "\n" ;

            } else { #inconsistence call
                print $TGT_union_conflict join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) . "\t" . join (", " , @$gt_all_in_one_locus) ."\t$reason_drop\n"  ;
            }
        }
        close $TGT_union_conflict ;
        close $TGT_union_vcf ;

        tabix_vcf($file_vcf_union_set) ;

        $file_vcfgz_foranalysis = "$file_vcf_union_set.gz" ;

    } elsif ($type_combine eq "intersection") {

        my $vcf_parser = Vcf_parser->new(file=> $file_vcfgz_foranalysis) ;
        my $file_output_isect_conflict = "$dir_results_output/multicaller_intersect_inconsistant_$jobid.txt" ;
        my $file_vcf_isect_set = "$dir_results_output/multicaller_insectset_$jobid.vcf" ;

        open (my $TGT_isect_conflict , ">$file_output_isect_conflict") ;
        open (my $TGT_isect_vcf , ">$file_vcf_isect_set") ;

        print $TGT_isect_vcf join("\t" , (qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/ , sort keys %$sample_list)) . "\n" ;

        while(my $var=$vcf_parser->next_var) {
            my $flag_inanalysis = 1 ;
            my $reason_drop = '' ;
            my $gt_all_in_one_locus = [] ;
            my $list_gt_by_sample = [] ;

            foreach my $sample (sort keys %$sample_list) {
                my $gt_got = {} ;
                foreach my $caller (sort keys %$caller_list) {
                    my $sample_tag = $sample . "_" . $caller ;
                    my $gt = $var->{sample_val}->{$sample_tag}->{GT} ;
                    $gt =~ s/\|/\// ;   # ignore phased, treat as an un-phased genotype
                    $gt_got->{$gt} ++ ;
                    push @$gt_all_in_one_locus, "$sample_tag:$gt" ;
                }

                my $list_of_got_gt = [keys %$gt_got] ;
                my $gt_type_count = scalar @$list_of_got_gt ;

                if ($gt_type_count == 1) {
                    $flag_inanalysis *= 1 ;
                    push @$list_gt_by_sample , $list_of_got_gt->[0] ;

                } else {
		    # Conflict genotype in 1+ callers of sample $sample  in insect
                     $flag_inanalysis *= 0 ;

                        $reason_drop .= "Conflict genotype in 1+ callers of sample $sample," ;
                }
            }
           if ($flag_inanalysis) {

                print $TGT_isect_vcf join("\t" , map {$var->{$_} } qw/CHROM POS ID REF ALT QUAL FILTER/ ) . "\t" . "\tGT\t" . join("\t",@$list_gt_by_sample) . "\n" ;

            } else {
                print $TGT_isect_conflict join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) . "\t" . join (", " , @$gt_all_in_one_locus) . "\t$reason_drop\n"  ;
            }
        }
        close $TGT_isect_conflict ;
        close $TGT_isect_vcf ;

        tabix_vcf($file_vcf_isect_set) ;


        $file_vcfgz_foranalysis = "$file_vcf_isect_set.gz" ;
    }
} # the end of multi_caller mode


# Check GT frequency
$log->write("Genotype frequency calculating start") ;
my $file_output_gtfreq = "$dir_results_output/aff_gt_freq_$jobid.txt" ;
my $file_vcfgz_gtfreq  = "$dir_results_output/aff_gt_freq_$jobid.vcf.gz" ;
my $file_output_gtfreq_unaff = "$dir_results_output/unaff_gt_freq_$jobid.txt" ;
my $file_vcfgz_gtfreq_unaff  = "$dir_results_output/unaff_gt_freq_$jobid.vcf.gz" ;

open (my $TGT_gtfreq ,">$file_output_gtfreq") ;
open (my $TGT_gtfreq_unaff ,">$file_output_gtfreq_unaff") ;

my $vcf_parser = Vcf_parser->new(file=>$file_vcfgz_foranalysis) ;
while (my $var = $vcf_parser->next_var) {
    my $affect_sample_num = scalar @$affect_samples ;
    my $unaffect_sample_num = scalar @$unaffect_samples ;

    if ($affect_sample_num) {
	my $gt_count = {} ;
	my $list_output_gt = [] ;
	my $list_output_freq = [] ;

	foreach my $sample (@$affect_samples) {
	    my $genotype = $var->{sample_val}->{$sample}->{GT} ;
	    $gt_count->{$genotype} ++ unless ($genotype eq '0/0'|| $genotype eq '0|0' || $genotype eq '0') ; # record , but skip ref-gt
	}
	my $gt_freq = {map {$_ => $gt_count->{$_} / $affect_sample_num } keys %$gt_count } ;
	my $alt_sum = 0 ;
	foreach my $gt (sort keys %$gt_freq) {
	    push @$list_output_gt   , $gt ;
	    push @$list_output_freq , $gt_freq->{$gt} ;
	    $alt_sum += $gt_freq->{$gt} ;
	}

	my $output_count = scalar @$list_output_gt ;
	if ($output_count) {
	    print $TGT_gtfreq join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) ;
	    print $TGT_gtfreq "\t" . join("," , @$list_output_gt) ;
	    print $TGT_gtfreq "\t" . join("," , @$list_output_freq) ;
	    print $TGT_gtfreq "\t$alt_sum" ;
	    print $TGT_gtfreq "\n" ;
	}
    }

    if ($unaffect_sample_num) {
	my $gt_count = {} ;
	my $list_output_gt_unaff = [] ;
	my $list_output_freq_unaff = [] ;

	foreach my $sample (@$unaffect_samples) {
            my $genotype = $var->{sample_val}->{$sample}->{GT} ;
            $gt_count->{$genotype} ++ unless ($genotype eq '0/0'|| $genotype eq '0|0' || $genotype eq '0') ; # record , but skip ref-gt
        }
	my $gt_freq = {map {$_ => $gt_count->{$_} / $unaffect_sample_num } keys %$gt_count } ;

	my $alt_sum =0 ;
	foreach my $gt (sort keys %$gt_freq) {
                push @$list_output_gt_unaff   , $gt ;
                push @$list_output_freq_unaff , $gt_freq->{$gt} ;
		$alt_sum += $gt_freq->{$gt} ;
	}

	my $output_count = scalar @$list_output_gt_unaff ;
	if ($output_count) {
	    print $TGT_gtfreq_unaff join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) ;
	    print $TGT_gtfreq_unaff "\t" . join("," , @$list_output_gt_unaff) ;
	    print $TGT_gtfreq_unaff "\t" . join("," , @$list_output_freq_unaff) ;
	    print $TGT_gtfreq_unaff "\t$alt_sum" ;
	    print $TGT_gtfreq_unaff "\n" ;
	}

    }
}

close $TGT_gtfreq ;
close $TGT_gtfreq_unaff ;

tabix_vcf($file_output_gtfreq,$log) ;
tabix_vcf($file_output_gtfreq_unaff,$log) ;

my $cmd_vannot_var_freq = "zcat $file_vcfgz_foranalysis| vcf-annotate " ;
$cmd_vannot_var_freq .= " -a $file_output_gtfreq.gz " ;
$cmd_vannot_var_freq .= " -c " . join("," , qw"CHROM POS REF ALT INFO/affected_vgt INFO/affected_vgfreq INFO/affected_altfreqsum") ;
$cmd_vannot_var_freq .= " -d key=INFO,ID=affected_vgt,Number=.,Type=String,Description='Variant genotype in affected samples' " ;
$cmd_vannot_var_freq .= " -d key=INFO,ID=affected_vgfreq,Number=.,Type=String,Description='Variant genotype frequency in affected samples' " ;
$cmd_vannot_var_freq .= " -d key=INFO,ID=affected_altfreqsum,Number=.,Type=Float,Description='Summary of affected-alt-frequency variant genotype in affected samples' " ;
$cmd_vannot_var_freq .= " 2> $dir_log/stderr_vannot_gtfreq_$jobid.log " ;
$cmd_vannot_var_freq .= " | bgzip -c > $file_vcfgz_gtfreq " ;

$log->andRun($cmd_vannot_var_freq) ;
tabix_vcf($file_vcfgz_gtfreq) ;

my $cmd_vannot_var_freq_unaff = "zcat $file_vcfgz_gtfreq| vcf-annotate " ;
$cmd_vannot_var_freq_unaff .= " -a $file_output_gtfreq_unaff.gz " ;
$cmd_vannot_var_freq_unaff .= " -c " . join("," , qw"CHROM POS REF ALT INFO/unaffected_vgt INFO/unaffected_vgfreq INFO/unaffected_altfreqsum") ;
$cmd_vannot_var_freq_unaff .= " -d key=INFO,ID=unaffected_vgt,Number=.,Type=String,Description='variant genotype in unaffected samples' " ;
$cmd_vannot_var_freq_unaff .= " -d key=INFO,ID=unaffected_vgfreq,Number=.,Type=String,Description='Variant genotype frequency in unaffected samples' " ;
$cmd_vannot_var_freq_unaff .= " -d key=INFO,ID=unaffected_altfreqsum,Number=.,Type=Float,Description='Summary of unaffected-alt-frequency variant genotype in affected samples' " ;
$cmd_vannot_var_freq_unaff .= " 2> $dir_log/stderr_vannot_gtfreq_$jobid.log " ;
$cmd_vannot_var_freq_unaff .= " | bgzip -c > $file_vcfgz_gtfreq_unaff " ;

$log->andRun($cmd_vannot_var_freq_unaff) ;
tabix_vcf($file_vcfgz_gtfreq_unaff) ;

my $collist_gtfreq = [qw/affected_vgt affected_vgfreq affected_altfreqsum unaffected_vgt unaffected_vgfreq unaffected_altfreqsum/] ;
my $type_list_gtfreq = [qw/text text float text text float/] ;
my $op_list_gtfreq = [qw/list list last list list last/] ;
my $gannot_list = {
    extract_list   => $collist_gtfreq ,
    type_list      => $type_list_gtfreq ,
    coladd_list    => [map {"$_\_$jobid"} @$collist_gtfreq] ,
    operation_list => $op_list_gtfreq ,
} ;

$log->write("Genotype frequency calculating finish") ;


my $vcf_input_for_next_step = $file_vcfgz_gtfreq_unaff ;

# if -c is set 
if ($file_cnv_list) {
    $log->write("CNV annotate start") ;

    my $file_vcfgz_cnv = "$dir_results_output/$jobid\_cnv.vcf.gz" ;
    my $cmd_cnv = "$dir_script/cnvkit_parser.pl " ;
    $cmd_cnv .= " -d $file_db " ;
    $cmd_cnv .= " -p $file_ped " ;
    $cmd_cnv .= " -c $file_cnv_list " ;
    $cmd_cnv .= " -i $vcf_input_for_next_step " ;
    $cmd_cnv .= " -o $file_vcfgz_cnv " ;
    $cmd_cnv .= " -j $jobid" ;

    $log->andRun($cmd_cnv) ;

    # prepare for gemini annotate
    my $collist_cnvkit = [qw/cnv_samples cnv_log2 cnv_fc_samples cnv_foldchange_log2/] ;
    push @{$gannot_list->{extract_list}}   , @$collist_cnvkit ;
    push @{$gannot_list->{type_list}}      , map {'text'} @$collist_cnvkit ;
    push @{$gannot_list->{coladd_list}}    , map {"$_\_$jobid"} @$collist_cnvkit ;
    push @{$gannot_list->{operation_list}} , map {'last'} @$collist_cnvkit ;
    
    $vcf_input_for_next_step = $file_vcfgz_cnv ;

    $log->write("CNV annotate finish") ;
}

# if -x is set
if ($file_xpr_list) {
    $log->write("expression profile parse and annotate start") ;
    my $script_xprparer = "$dir_script/xprprofile_parser.pl" ;
    my $file_vcfgz_xpr = "$file_prefix\_xpr.vcf.gz" ;

    my $cmd_xpr = "$script_xprparer " ;
    $cmd_xpr .= " -p $file_ped " ;
    $cmd_xpr .= " -x $file_xpr_list " ;
#    $cmd_xpr .= " -i $file_vcfgz_foranalysis " ;
    $cmd_xpr .= " -i $vcf_input_for_next_step " ;
    $cmd_xpr .= " -o $file_vcfgz_xpr " ;
    $cmd_xpr .= " -j $jobid " ;
    $cmd_xpr .= " -l $dir_log " ;
    $log->andRun($cmd_xpr) ;

    # prepare for gemini annotate
    my $collist_xpr = [qw/xpr_samples xpr_readcount xpr_tpm xpr_fpkm xpr_exonreadcount xpr_exonid xpr_fc_samples xpr_foldchange_readcount xpr_foldchange_tpm xpr_foldchange_fpkm xpr_foldchange_exonreadcount/] ;

    push @{$gannot_list->{extract_list}} , @$collist_xpr ;
    push @{$gannot_list->{type_list}}      , map {'text'} @$collist_xpr ;
    push @{$gannot_list->{coladd_list}}    , (qw/xpr_samples xpr_readcount xpr_tpm xpr_fpkm xpr_exonreadcount xpr_exonid xpr_fc_samples/ ,  "xpr_foldchange_readcount_$jobid", "xpr_foldchange_tpm_$jobid", "xpr_foldchange_fpkm_$jobid", "xpr_foldchange_exonreadcount_$jobid" ) ;
    push @{$gannot_list->{operation_list}} , map {'last'} @$collist_xpr ;

    $vcf_input_for_next_step = $file_vcfgz_xpr ;

    $log->write("expression profile parse and annotate finish") ;

}

# -m = paired or family
my $col_list_for_wf_union_set = [] ;
if ($mode_workflow eq 'family') {
    my $file_vcf_for_workflow = "$dir_work/for_genetic_workflow_$jobid.vcf" ;
    my $cmd_extract_for_workflow = "zcat $vcf_input_for_next_step > $file_vcf_for_workflow" ;
    $log->andRun($cmd_extract_for_workflow) ;

    my $file_stderr_wf_ar = "$dir_log/stderr_runworkflow_AR_$jobid.log" ;
    my $file_stderr_wf_ch = "$dir_log/stderr_runworkflow_CH_$jobid.log" ;
    my $file_stderr_wf_dr = "$dir_log/stderr_runworkflow_DR_$jobid.log" ;
    my $file_stderr_wf_th = "$dir_log/stderr_runworkflow_TH_$jobid.log" ;
    my $file_stderr_wf_xl = "$dir_log/stderr_runworkflow_XL_$jobid.log" ;

    my $file_gz_output_ar = "$dir_results_output/workflow_AR_output_$jobid.txt.gz" ;
    my $file_gz_output_ch = "$dir_results_output/workflow_CH_output_$jobid.txt.gz" ;
    my $file_gz_output_dr = "$dir_results_output/workflow_DR_output_$jobid.txt.gz" ;
    my $file_gz_output_th = "$dir_results_output/workflow_TH_output_$jobid.txt.gz" ;
    my $file_gz_output_xl = "$dir_results_output/workflow_XL_output_$jobid.txt.gz" ;

    my $col_ar = "is_AR_$jobid" ;
    my $col_ch = "is_CH_$jobid" ;
    my $col_dr = "is_DR_$jobid" ;
    my $col_th = "is_TH_$jobid" ;
    my $col_xl = "is_XL_$jobid" ;

    my $file_vcfgz_ar_output = "$dir_results_output/AutosomalRecessive_$jobid.vcf.gz" ;
    my $file_vcfgz_ch_output = "$dir_results_output/CompoundHet_wAr_$jobid.vcf.gz" ;
    my $file_vcfgz_dr_output = "$dir_results_output/DenovoRecessive_wArCh_$jobid.vcf.gz" ;
    my $file_vcfgz_th_output = "$dir_results_output/TwoHits_wArChDr_$jobid.vcf.gz" ;
    my $file_vcfgz_xl_output = "$dir_results_output/XLinked_wArChDrTh_$jobid.vcf.gz" ;

    push @$col_list_for_wf_union_set, map { $_.$jobid } qw/is_AR_ is_CH_ is_DR_ is_TH_ is_XL_/ ;

    # Autosomal-recessive
    my $exe_wf_ar  = "$dir_workflows/autosomal-recessive/Autosomal-recessive.py" ;
    my $exe_wf2_ar = "$dir_workflows/autosomal-recessive/Compare.py" ;
    running_workflow($exe_wf_ar , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_ar , $file_gz_output_ar , $exe_wf2_ar) ;
    wf_vannot($file_gz_output_ar , $col_ar , $vcf_input_for_next_step , $file_vcfgz_ar_output ) ;

    # Compound-het
    my $exe_wf_ch  = "$dir_workflows/compound-het/Compound-het.py" ;
    my $exe_wf2_ch = "$dir_workflows/compound-het/Compare.py" ;
    running_workflow($exe_wf_ch , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_ch , $file_gz_output_ch , $exe_wf2_ch , $file_db) ;
    wf_vannot($file_gz_output_ch , $col_ch , $file_vcfgz_ar_output , $file_vcfgz_ch_output) ;

    # De Novo-recessive
    my $exe_wf_dr  =  "$dir_workflows/denovo-recessive/Denovo-recessive.py" ;
    my $exe_wf2_dr =  "$dir_workflows/denovo-recessive/Compare.py" ;
    running_workflow($exe_wf_dr , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_dr , $file_gz_output_dr,$exe_wf2_dr) ;
    wf_vannot($file_gz_output_dr, $col_dr , $file_vcfgz_ch_output , $file_vcfgz_dr_output) ;

    # Two-hit
    my $exe_wf_th  = "$dir_workflows/two-hits/Two-hits.py" ;
    my $exe_wf2_th = "$dir_workflows/two-hits/Compare.py" ;
    running_workflow($exe_wf_th , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_th , $file_gz_output_th , $exe_wf2_th) ;
    wf_vannot($file_gz_output_th , $col_th , $file_vcfgz_dr_output , $file_vcfgz_th_output) ;

    # X-linked
    my $exe_wf_xl  = "$dir_workflows/x-linked/X-linked.py" ;
    my $exe_wf2_xl = "$dir_workflows/x-linked/Compare.py" ;
    running_workflow($exe_wf_xl , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_xl , $file_gz_output_xl , $exe_wf2_xl) ;
    wf_vannot($file_gz_output_xl , $col_xl , $file_vcfgz_th_output , $file_vcfgz_xl_output) ;

    $vcf_input_for_next_step = $file_vcfgz_xl_output ;

    my $in_analysis_file_list = {
	$file_gz_output_ar => 5,
	$file_gz_output_ch => 5,
	$file_gz_output_dr => 5,
	$file_gz_output_th => 5,
	$file_gz_output_xl => 5,
    } ;

    my $file_gz_inanalysis = "$dir_results_output/inanalysis_$jobid.txt.gz" ;
    my $file_vcfgz_inanalysis = "$dir_results_output/inanalysis_$jobid.vcf.gz" ;

    get_union($in_analysis_file_list , $file_gz_inanalysis ) ;
    tabix_vcf($file_gz_inanalysis) ;

    my $cmd_vannot_inanalysis = "zcat $vcf_input_for_next_step | vcf-annotate " ;
    $cmd_vannot_inanalysis .= " -a $file_gz_inanalysis " ;
    $cmd_vannot_inanalysis .= " -c " . join("," , ('CHROM', 'POS','REF','ALT',"INFO/in_analysis_$jobid")) ;
    $cmd_vannot_inanalysis .= " -d key=INFO,ID=in_analysis_$jobid,Number=1,Type=Integer,Description='in_analysis_$jobid' " ;
    $cmd_vannot_inanalysis .= "2> $dir_log/stderr_vannot_inanalysis_$jobid.log " ;
    $cmd_vannot_inanalysis .= " |bgzip -c > $file_vcfgz_inanalysis " ;
    $log->andRun($cmd_vannot_inanalysis) ;
    
    $vcf_input_for_next_step = $file_vcfgz_inanalysis ;

    push @{$gannot_list->{extract_list}} , ($col_ar,$col_ch,$col_dr,$col_th,$col_xl,"in_analysis_$jobid") ;
    push @{$gannot_list->{type_list}} , qw/integer integer integer integer integer integer/ ;
    push @{$gannot_list->{coladd_list}} , ($col_ar,$col_ch,$col_dr,$col_th,$col_xl,"in_analysis_$jobid") ;
    push @{$gannot_list->{operation_list}} , qw/last last last last last last/ ;

} elsif ($mode_workflow eq 'paired') {

    # LOH detect
    $log->write("LOH detect start") ;
    my $script_lohdetect = "$dir_script/loh_detector.pl" ;
    my $file_vcfgz_isLOH = "$file_prefix\_loh.vcf.gz" ;

    my $cmd_loh = "$script_lohdetect" ;
    $cmd_loh .= " -j $jobid " ;
#    $cmd_loh .= " -i $file_vcfgz_foranalysis " ;
    $cmd_loh .= " -i $vcf_input_for_next_step " ;
    $cmd_loh .= " -o $file_vcfgz_isLOH " ;
    $cmd_loh .= " -p $file_ped " ;
    $cmd_loh .= " -d $file_db " ;
    $cmd_loh .= " -n $threads_enable " ;

    $log->andRun($cmd_loh) ;
    $log->loadlog("$dir_log/log_lohdetect_$jobid.log") ;
    $log->write("LOH detect finish") ;
    push @$col_list_for_wf_union_set, "is_LOH_$jobid" ;

    # Somatic detect
    $log->write("Somatic detect start") ;
    my $script_somdetect = "$dir_script/somatic_detector.pl" ;
    my $file_output_isSomatic = "$dir_results_output/is_Somatic_$jobid.txt" ;
    my $file_vcfgz_isSomatic = "$file_prefix\_somatic.vcf.gz" ;

    my $cmd_somatic = "$script_somdetect " ;
#    $cmd_somatic .= " -i $file_vcfgz_foranalysis " ;
    $cmd_somatic .= " -i $file_vcfgz_isLOH " ;
    $cmd_somatic .= " -o $file_output_isSomatic " ;
    $cmd_somatic .= " -v $file_vcfgz_isSomatic " ;
    $cmd_somatic .= " -p $file_ped " ;
    $cmd_somatic .= " -j $jobid " ;
    $cmd_somatic .= " -d $file_db " ;
    $cmd_somatic .= " -n $threads_enable " ;

    $log->andRun($cmd_somatic) ;

    tabix_vcf("$file_output_isSomatic") ;

    $log->write("Somatic detect finish") ;
    push @$col_list_for_wf_union_set, "is_Somatic_$jobid" ;
    

    my $collist_pairworkflow = [qw/is_Somatic Somatic_gt is_LOH in_analysis/] ;
    my $typelist_pairworkflow = [qw/integer text integer integer/] ;
    push @{$gannot_list->{extract_list}}   , @$collist_pairworkflow ;
    push @{$gannot_list->{type_list}}      , @$typelist_pairworkflow ;
    push @{$gannot_list->{coladd_list}}    , map {"$_\_$jobid"} @$collist_pairworkflow ;
    push @{$gannot_list->{operation_list}} , map {'last'} @$collist_pairworkflow ;

    $vcf_input_for_next_step = $file_vcfgz_isSomatic ;
    
} else {
    # none
}

store $gannot_list , $file_hashref_gannot_tmp ;

my $cmd_gannot = "$dir_script/gannot.pl -j $jobid -d $file_db -g $file_hashref_gannot_tmp -v $vcf_input_for_next_step " ;
$log->andRun($cmd_gannot) ;

close $CONFIG ;
$log->write("VSanlz finish") ;

sub wf_vannot {
    my $file_annotation = shift ;
    my $tag = shift ;
    my $file_vcf_input = shift ;
    my $file_vcf_output = shift ;

    my $cmd_vannot = "vcf-annotate " ;
    $cmd_vannot .= " -a $file_annotation " ;
    $cmd_vannot .= " -c " . join("," , ('CHROM', 'POS','REF','ALT',"INFO/$tag")) ;
    $cmd_vannot .= " -d key=INFO,ID=$tag,Number=1,Type=Integer,Description='genetic model $tag' " ;
    $cmd_vannot .= "$file_vcf_input 2> $dir_log/stderr_vannot_workflow_$tag.log " ;
    $cmd_vannot .= " |bgzip -c > $file_vcf_output " ;

    $log->andRun($cmd_vannot);

    tabix_vcf($file_vcf_output) ;
}

sub running_workflow {
    my $script = shift ;
    my $file_vcf = shift ;
    my $file_ped = shift ;
    my $file_stderr = shift ;
    my $file_gz_output_wf = shift ;

    my $script_compare = shift ;
    my $db = shift || 0 ;

    my $output_tmp = "$file_gz_output_wf.tmp" ;

    # workflow step 1 of genetic model
    my $cmd_runworkflow = "$script " ;
    $cmd_runworkflow .= " -v $file_vcf " ;
    $cmd_runworkflow .= " -p $file_ped " ;
    $cmd_runworkflow .= " -d $db " if ($db) ;
    $cmd_runworkflow .= " 2> $file_stderr " ;
    $cmd_runworkflow .= " > $output_tmp" ;

    $log->andRun($cmd_runworkflow) ;

    if (-e $script_compare) {
	# workflow step 2 of genetic model
        my $cmd_step2 = "$script_compare " ;
        $cmd_step2 .= " -v $file_vcf " ;
        $cmd_step2 .= " -p $file_ped " ;
        $cmd_step2 .= " -c $output_tmp " ;
        $cmd_step2 .= " | bgzip -c > $file_gz_output_wf" ;
        $log->andRun($cmd_step2) ;

    } else { # without Compare.py, run step 1 ONLY
        my $cmd_bgzip_output_tmp = "bgzip -@ $threads_enable -c $output_tmp > $file_gz_output_wf" ;
        $log->andRun($cmd_bgzip_output_tmp) ;
    }


    tabix_vcf($file_gz_output_wf) ;
}


sub get_union {
    my $files = shift ;
    my $file_union = shift ;

    my $union_set = {} ;
    foreach my $file (keys %$files) {

	my $col1 = $files->{$file} ;

	open (my $SRC,"zcat $file|") ;
	while(<$SRC>) {
	    chomp ;
	    my @data = split /\t/ ;
	    my $col = $col1 - 1 ;

	    if ($data[$col]) {
		my $tag = join("\t" , ($data[0],$data[1],$data[2],$data[3])) ;

		$union_set->{$tag} = 1 ;

	    }
	}
	close $SRC ;
    }

    open (my $TGT,">$file_union.tmp") ;

    foreach my $tag (keys %$union_set) {
	print $TGT "$tag\t1\n" ;
    }

    close $TGT ;

    `vcf-sort $file_union.tmp |bgzip -c > $file_union` ;
}
