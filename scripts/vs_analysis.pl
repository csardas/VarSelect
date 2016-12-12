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
my $file_cnv_list = $opts->{c} ;
my $file_xpr_list = $opts->{x} ;

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
open (my $SRC_config , "$file_config_db" ) ;
while (<$SRC_config>) {
    chomp ;
    my($idx,$val) = split /\t/ ;
    $config_from_db->{$idx} = $val ;
}
close ;

open (my $CONFIG , ">$file_config") ;
print $CONFIG "ID\t$jobid\n" ;
print $CONFIG "DB\t$file_db\n" ;

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

# if multicaller mode, get combined genotype between callers,  generate variant list with inconsistence genotype

# Dump VCF from gemini db 
my $file_vcfgz_foranalysis = "$fpath_db/$fname_db.vcf.gz" ; # default setting

if (-e $file_vcfgz_foranalysis) {
    tabix_vcf($file_vcfgz_foranalysis,$log) ; # unless (-e "$file_vcfgz_foranalysis.tbi") ;

} else {
    die("$file_vcfgz_foranalysis not found! please check your db create directory and put the $fname_db.vcf.gz to the current location of your db!\n") ;
}

if ($flag_multicallers) {
    if ($type_combine eq "union") {

        my $vcf_parser = Vcf_parser->new(file=> $file_vcfgz_foranalysis) ;
        my $file_output_union_temp = "$dir_log/union_tmp_$jobid.txt" ;

        my $file_output_union_conflict = "$dir_log/union_conflict_$jobid.txt" ;
        my $file_vcf_union_set = "$dir_log/$file_prefix\_multicaller_unionset_$jobid.vcf";

        open (my $TGT_union_tmp , ">$file_output_union_temp") ;
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
                # output bed format for Gemini annotate
                print $TGT_union_tmp join("\t" , map {$var->{$_} } qw/CHROM POS POS/ ) . "\t1\n" ;    # bed like format

                # output vcf with sample only (no caller now) for workflow analysis
                 print $TGT_union_vcf join("\t" , map {$var->{$_} } qw/CHROM POS ID REF ALT QUAL FILTER/ ) . "\tin_analysis_$jobid;" . "\tGT\t" . join("\t",@$list_gt_by_sample) . "\n" ;
            } else { #inconsistence call
                print $TGT_union_conflict join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) . "\t" . join (", " , @$gt_all_in_one_locus) ."\t$reason_drop\n"  ;
            }
        }
       close $TGT_union_tmp ;
        close $TGT_union_conflict ;
        close $TGT_union_vcf ;

        tabix_vcf($file_vcf_union_set) ;

        bgzip($file_output_union_temp, "$file_output_union_temp.gz") ;

        my $cmd_tabix_uniontmp = "tabix -f -s 1 -b 2 -e 2 -c # $file_output_union_temp.gz" ;
        $log->write($cmd_tabix_uniontmp) ;
        `$cmd_tabix_uniontmp` ;

        my $cmd_gannot = "gemini annotate " ;
        $cmd_gannot .= " -a boolean " ;
        $cmd_gannot .= " -f $file_output_union_temp.gz " ;
        $cmd_gannot .= " -c in_analysis_$jobid " ;
        $cmd_gannot .= " $file_db 2> $dir_log/stderr_analysis_gannot_union_$jobid.log" ;

        $log->write($cmd_gannot) ;
        `$cmd_gannot` ;

        $file_vcfgz_foranalysis = "$file_vcf_union_set.gz" ;

    } elsif ($type_combine eq "intersection") {

        my $vcf_parser = Vcf_parser->new(file=> $file_vcfgz_foranalysis) ;
        my $file_output_isect_temp = "$dir_log/isect_tmp_$jobid.txt" ;
        my $file_output_isect_conflict = "$dir_log/isect_conflict_$jobid.txt" ;
        my $file_vcf_isect_set = "$dir_log/$file_prefix\_multicaller_insectset_$jobid.vcf" ;

        open (my $TGT_isect_tmp , ">$file_output_isect_temp") ;
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
                print $TGT_isect_tmp join("\t" , map {$var->{$_} } qw/CHROM POS POS/ ) ;         #bed format
                print $TGT_isect_tmp "\t1\n" ;

                print $TGT_isect_vcf join("\t" , map {$var->{$_} } qw/CHROM POS ID REF ALT QUAL FILTER/ ) . "\tin_analysis_$jobid;" . "\tGT\t" . join("\t",@$list_gt_by_sample) . "\n" ;

            } else {
                print $TGT_isect_conflict join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) . "\t" . join (", " , @$gt_all_in_one_locus) . "\t$reason_drop\n"  ;
            }
        }
        close $TGT_isect_tmp ;
        close $TGT_isect_conflict ;
        close $TGT_isect_vcf ;

        tabix_vcf($file_vcf_isect_set) ;

        bgzip ($file_output_isect_temp , "$file_output_isect_temp.gz") ;
        my $cmd_tabix_isect_tmp = "tabix -f -s 1 -b 2 -e 2 -c # $file_output_isect_temp.gz " ;
        $log->andRun($cmd_tabix_isect_tmp) ;

        my $cmd_gannot = "gemini annotate " ;
        $cmd_gannot .= " -a boolean " ;
        $cmd_gannot .= " -f $file_output_isect_temp.gz " ;
        $cmd_gannot .= " -c in_analysis_$jobid " ;
        $cmd_gannot .= " $file_db 2> $dir_log/stderr_analysis_gannot_isect_$jobid.log" ;
        $log->andRun($cmd_gannot) ;

        $file_vcfgz_foranalysis = "$file_vcf_isect_set.gz" ;
    }


}

# if -c is set 
if ($file_cnv_list) {
    $log->write("CNV annotate start") ;

    my $file_vcfgz_cnv = "$dir_results_output/$jobid\_cnv.vcf.gz" ;
    my $cmd_cnv = "$dir_script/cnvkit_parser.pl -d $file_db -p $file_ped -c $file_cnv_list -i $file_vcfgz_foranalysis -o $file_vcfgz_cnv -j $jobid" ;
    $log->andRun($cmd_cnv) ;
#    $log->loadlog("$dir_log/log_cnvkit_parser_$jobid.log") ;
    $log->write("CNV annotate finish") ;
}

# if -x is set
if ($file_xpr_list) {
    $log->write("dexseq parse and annotate start") ;
    my $script_dexseqparer = "$dir_script/dexseq_parser.pl" ;
    my $file_vcfgz_xpr = "$file_prefix\_xpr.vcf.gz" ;

    my $cmd_xpr = "$script_dexseqparer -p $file_ped -x $file_xpr_list -i $file_vcfgz_foranalysis -o $file_vcfgz_xpr -j $jobid -l $dir_log " ;
    $log->andRun($cmd_xpr) ;

    my $cmd_gannot_xpr = "gemini annotate -f $file_vcfgz_xpr " ;
    $cmd_gannot_xpr .= " -a extract " ;
    $cmd_gannot_xpr .= " -e " . join("," , qw/xpr_samples xpr_readcount xpr_tpm xpr_fpkm xpr_exonreadcount xpr_exonid xpr_fc_samples xpr_foldchange_readcount xpr_foldchange_tpm xpr_foldchange_fpkm xpr_foldchange_exonreadcount/) ;
    $cmd_gannot_xpr .= " -t " . join("," , qw/text text text text text text text text text text text/) ;
    $cmd_gannot_xpr .= " -c " . join("," , ( qw/xpr_samples xpr_readcount xpr_tpm xpr_fpkm xpr_exonreadcount xpr_exonid xpr_fc_samples/ , "xpr_foldchange_readcount_$jobid", "xpr_foldchange_tpm_$jobid", "xpr_foldchange_fpkm_$jobid", "xpr_foldchange_exonreadcount_$jobid") ) ;
    $cmd_gannot_xpr .= " -o " . join("," , qw/first first first first first first first first first first first/) ;
    $cmd_gannot_xpr .= " $file_db 2> $dir_log/stderr.gannot_xpr.jobid.log" ;
    $log->andRun($cmd_gannot_xpr) ;

    $log->write("dexseq parse and annotate finish") ;

}

# -m = paired or family
my $col_list_for_wf_union_set = [] ;
if ($mode_workflow eq 'family') {
    my $file_vcf_for_workflow = "$dir_work/for_genetic_workflow_$jobid.vcf" ;
    my $cmd_extract_for_workflow = "zcat $file_vcfgz_foranalysis > $file_vcf_for_workflow" ;
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

    my $file_stderr_gannot_ar = "$dir_log/stderr_gannot_AR_$jobid.log" ;
    my $file_stderr_gannot_ch = "$dir_log/stderr_gannot_CH_$jobid.log" ;
    my $file_stderr_gannot_dr = "$dir_log/stderr_gannot_DR_$jobid.log" ;
    my $file_stderr_gannot_th = "$dir_log/stderr_gannot_TH_$jobid.log" ;
    my $file_stderr_gannot_xl = "$dir_log/stderr_gannot_XL_$jobid.log" ;


    push @$col_list_for_wf_union_set, map { $_.$jobid } qw/is_AR_ is_CH_ is_DR_ is_TH_ is_XL_/ ;

    # Autosomal-recessive
    my $exe_wf_ar = "$dir_workflows/autosomal-recessive/Autosomal-recessive.py" ;
    my $exe_wf2_ar = "$dir_workflows/autosomal-recessive/Compare.py" ;
    running_workflow($exe_wf_ar , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_ar , $file_gz_output_ar , $exe_wf2_ar) ;
    wf_vgannot($file_gz_output_ar,$col_ar,$file_vcfgz_foranalysis , "$dir_results_output/AutosomalRecessive_$jobid.vcf.gz" , $file_db ) ;

    # Compound-het
    my $exe_wf_ch = "$dir_workflows/compound-het/Compound-het.py" ;
    my $exe_wf2_ch = "$dir_workflows/compound-het/Compare.py" ;
    running_workflow($exe_wf_ch , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_ch , $file_gz_output_ch , $exe_wf2_ch , $file_db) ;
    wf_vgannot($file_gz_output_ch,$col_ch,$file_vcfgz_foranalysis , "$dir_results_output/CompoundHet_$jobid.vcf.gz" , $file_db ) ;

    # De Novo-recessive
    my $exe_wf_dr =  "$dir_workflows/denovo-recessive/Denovo-recessive.py" ;
    my $exe_wf2_dr =  "$dir_workflows/denovo-recessive/Compare.py" ;
    running_workflow($exe_wf_dr , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_dr , $file_gz_output_dr,$exe_wf2_dr) ;
    wf_vgannot($file_gz_output_dr,$col_dr,$file_vcfgz_foranalysis , "$dir_results_output/DenovoRecessive_$jobid.vcf.gz" , $file_db ) ;

    # Two-hit
    my $exe_wf_th = "$dir_workflows/two-hits/Two-hits.py" ;
    my $exe_wf2_th = "$dir_workflows/two-hits/Compare.py" ;
    running_workflow($exe_wf_th , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_th , $file_gz_output_th , $exe_wf2_th) ;
    wf_vgannot($file_gz_output_th,$col_th,$file_vcfgz_foranalysis , "$dir_results_output/TwoHits_$jobid.vcf.gz" , $file_db ) ;

    # X-linked
    my $exe_wf_xl = "$dir_workflows/x-linked/X-linked.py" ;
    my $exe_wf2_xl = "$dir_workflows/x-linked/Compare.py" ;
    running_workflow($exe_wf_xl , $file_vcf_for_workflow , $file_ped , $file_stderr_wf_xl , $file_gz_output_xl , $exe_wf2_xl) ;
    wf_vgannot($file_gz_output_xl,$col_xl,$file_vcfgz_foranalysis , "$dir_results_output/XLinked_$jobid.vcf.gz" , $file_db ) ;


} elsif ($mode_workflow eq 'paired') {

    # LOH detect
    $log->write("LOH detect start") ;
    my $script_lohdetect = "$dir_script/loh_detector.pl" ;
    my $file_vcfgz_isLOH = "$file_prefix\_loh.vcf.gz" ;

    my $cmd_loh = "$script_lohdetect" ;
    $cmd_loh .= " -j $jobid " ;
    $cmd_loh .= " -i $file_vcfgz_foranalysis " ;
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
    $cmd_somatic .= " -i $file_vcfgz_foranalysis " ;
    $cmd_somatic .= " -o $file_output_isSomatic " ;
    $cmd_somatic .= " -p $file_ped " ;
    $cmd_somatic .= " -j $jobid " ;
    $cmd_somatic .= " -d $file_db " ;
    $cmd_somatic .= " -n $threads_enable " ;

    $log->andRun($cmd_somatic) ;

    tabix_vcf("$file_output_isSomatic") ;

    $log->write("Somatic detect finish") ;
    push @$col_list_for_wf_union_set, "is_Somatic_$jobid" ;
    
} else {
    # should not be here
}

# get union set from all paired workflows
my $file_final_union_set = "$dir_results_output/final_union_set.txt" ;
my $gemini_final_union = "gemini query -q 'select chrom,start+1,ref,alt from variants where " ;

$gemini_final_union .= join(" OR " , map { "$_ = '1'" }  @$col_list_for_wf_union_set ) ;
$gemini_final_union .= "' $file_db |vcf-sort > $file_final_union_set" ;
$log->andRun($gemini_final_union) ;

bgzip($file_final_union_set , "$file_final_union_set.gz" ) ;
tabix_vcf("$file_final_union_set.gz" ) ;

# gemini annotate union set
my $cmd_gannot_final = "gemini annotate " ;
$cmd_gannot_final .= " -f $file_final_union_set.gz " ;
$cmd_gannot_final .= " -a boolean " ;
$cmd_gannot_final .= " -c in_analysis_$jobid " ;
$cmd_gannot_final .= "$file_db 2> $dir_log/stderr_gannot_final_$jobid.log" ;

$log->andRun($cmd_gannot_final) ;



close $CONFIG ;
$log->write("VSanlz finish") ;

sub wf_vgannot {
    my $file_annotation = shift ;
    my $tag = shift ;
    my $file_vcf_input = shift ;
    my $file_vcf_output = shift ;
    my $file_db = shift ;

    my $cmd_vannot = "vcf-annotate " ;
    $cmd_vannot .= " -a $file_annotation " ;
    $cmd_vannot .= " -c " . join("," , ('CHROM', 'POS','REF','ALT',"INFO/$tag")) ;
    $cmd_vannot .= " -d key=INFO,ID=$tag,Number=1,Type=Integer,Description='genetic model $tag' " ;
    $cmd_vannot .= "$file_vcf_input 2> $dir_log/stderr_vannot_workflow_$tag.log " ;
    $cmd_vannot .= " |bgzip -c > $file_vcf_output " ;

    $log->andRun($cmd_vannot);

    tabix_vcf($file_vcf_output) ;

    my $cmd_gannot = "gemini annotate " ;
    $cmd_gannot .= " -f $file_vcf_output " ;
    $cmd_gannot .= " -a extract " ;
    $cmd_gannot .= " -c $tag " ;
    $cmd_gannot .= " -e $tag " ;
    $cmd_gannot .= " -t integer " ;
    $cmd_gannot .= " -o first " ;
    $cmd_gannot .= " 2> $dir_log/stderr_wf_$tag.log " ;
    $cmd_gannot .= " $file_db " ;

    $log->andRun($cmd_gannot) ;
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

    my $cmd_runworkflow = "$script " ;
    $cmd_runworkflow .= " -v $file_vcf " ;
    $cmd_runworkflow .= " -p $file_ped " ;
    $cmd_runworkflow .= " -d $db " if ($db) ;
    $cmd_runworkflow .= " 2> $file_stderr " ;
    $cmd_runworkflow .= " > $output_tmp" ;

    $log->andRun($cmd_runworkflow) ;

    if (-e $script_compare) {
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

