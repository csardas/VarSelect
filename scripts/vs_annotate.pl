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

# default global setting
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

# options handle
my $opts = {} ;
getopts("v:p:j:n:kd:",$opts) ;

my $file_vcf_list	= $opts->{v} ;
my $file_ped		= $opts->{p} ;
my $jobid		= $opts->{j} ;
my $flag_multicallers   = $opts->{k} || 0 ;
my $file_db		= $opts->{d} ;
my $threads_enable      = $opts->{n} ;
my $threads_active = $threads_enable ;

# jobs specific path
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;
my $file_config = "$dir_results_output/varselect.config" ;
my $file_log = "$dir_log/log_vsannotate_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("VSannot start") ;
$log->write("Jobid: $jobid") ;
$log->write(join (" " , map {"-$_ $opts->{$_}"} keys %$opts ) ) ;

open (my $CONFIG , ">$file_config") ;
print $CONFIG "ID\t$jobid\n" ;

my $file_hashref_vannot_tmp = "$dir_work/vannot_tmp.hashref" ;

my ($fname_db,$fpath_db,$fext_db) = fileparse($file_db,qw/.db/) ;
my $file_geminidb = $file_db ;
my $file_db_prefix = "$fpath_db/$fname_db" ;
my $file_vcfgz_foranalysis = "$file_db_prefix.vcf.gz" ;
my $file_config_db = "$file_db_prefix.config" ;

my $file_prefix = "$dir_results_output/$fname_db" ;
my ($file_prefix0 , $file_path , $file_ext) = fileparse($file_vcf_list , qw/.txt .csv/) ;

# PreProcess 01: Load vcf files of samples and get sample/caller list
#=======================================
my $sample2file = {} ;
my $caller_list = {} ;
my $sample_list = {} ;
$log->write ("Load and check VCF file list start") ;
open (my $SRC_filelist , "$file_vcf_list") ;
while (my $line = <$SRC_filelist>) {
    chomp $line ;

    my ($sample, $caller, $file_vcf) ;
    my $sample_tag ;

    if ($flag_multicallers) {	# multi-caller mode
        ($sample, $caller, $file_vcf) = split (/\,/ , $line) ;

        die "multi-caller mode\n Please provide \"Sample, Caller, FileLocation\" foreach input file." unless (defined $file_vcf) ;
        $caller_list->{$caller} ++ ;
        $sample_tag = "$sample\_$caller" ;

    } else {			# single-caller mode
        ($sample, $file_vcf) = split (/\,/ , $line) ;
        $sample_tag = $sample ;
    }

    die ("VCF file $file_vcf is not exists.\n") unless (-e $file_vcf) ;
    $sample_list->{$sample} ++ ;

    my ($filevcf_fn , $filevcf_path, $filevcf_ext) = fileparse($file_vcf , qw/.vcf.gz .vcf/) ;

    # make sure every input vcf is bgziped and tabixed
    my $file_vcfbgz_input = "$dir_work/$filevcf_fn.vcf.gz" ;
    my $cmd_force_bgzipped_input_vcf = "$dir_script/drop_multiline.pl -j $jobid -i $file_vcf |bgzip -@ $threads_enable -c > $file_vcfbgz_input " ;
    $log->andRun($cmd_force_bgzipped_input_vcf) ;
    tabix_vcf($file_vcfbgz_input) ;

    tabix_vcf($file_vcfbgz_input) ;

    push @{$sample2file->{$sample_tag}} , $file_vcfbgz_input ;
}
close $SRC_filelist ;
$log->write ("Load and check VCF file list finish") ;

print $CONFIG "multicaller\t$flag_multicallers\n" ;
print $CONFIG "caller_list\t" . join("," , sort keys %$caller_list) . "\n" if ($flag_multicallers) ;
print $CONFIG "sample_list\t" . join("," , sort keys %$sample_list) . "\n" ;

# PreProcess 02: Load PED file, output Sample sex, and  get sample list
#=======================================
my $file_sample_sex = $file_prefix . "_sample_sex.txt" ;
my $file_sample_sex_mcaller = $file_prefix . "_sample_sex_mcaller.txt" ;
my $list_samples = {} ;

my ($ped_fn,$ped_path,$ped_ext) = fileparse($file_ped,qw/.ped/) ;
my $file_ped_mcaller = $ped_path . $ped_fn . "_mcaller.ped" ;
open (my $TGT_pedmcaller,">$file_ped_mcaller") if ($flag_multicallers) ;

$log->write ("Load ped file start") ;
open (my $TGT_samplesex , ">$file_sample_sex") ;
open (my $TGT_samplesex_mcaller , ">$file_sample_sex_mcaller") if ($flag_multicallers) ;

open (my $SRC_ped , "$file_ped") ;
while (<$SRC_ped>) {
    next if /^#/ ;  # skip line start with #
    chomp ;
    my ($family,$sample,$fatherid,$motherid,$sexid,$affectid,$other) = split /\t/ ;

    my $sex = '' ;
    if ($sexid == 1) {
        $sex = 'M' ;
    } elsif ($sexid == 2) {
        $sex = 'F' ;
    } else {
        die "PED file $file_ped with wrong sex id $sexid\n\n" ;
    }

    print $TGT_samplesex "$sample $sex\n" ;
    if ($flag_multicallers) {
	# alternate ped file
        foreach my $caller (sort keys %$caller_list) {
            my $sample_tag = $sample . "_" . $caller ;
            die "$sample_tag doesn't have related vcf file, please check." unless (exists $sample2file->{$sample_tag}) ;
            print $TGT_pedmcaller join("\t" , ($family,$sample_tag,($fatherid)?"$fatherid\_$caller":$fatherid , ($motherid)?"$motherid\_$caller":$motherid , $sex, $affectid, $other)) ;
            print $TGT_pedmcaller "\n" ;
	    # alternate sex file
	    print $TGT_samplesex_mcaller "$sample_tag $sex\n" ;
            $list_samples->{$sample_tag} ++ ;

            # reheader by bcftools sample => sample_caller
	    my $file_new_header = "$dir_work/newheader_$sample_tag.txt" ;
            open (my $TGT_nh,">$file_new_header") ;
            print $TGT_nh "$sample $sample_tag\n" ;
            close $TGT_nh ;

            foreach my $old_vcfgz (@{$sample2file->{$sample_tag}}) {
                my $file_new_vcf = "$dir_work/$sample_tag\_nh.vcf" ;
                my $cmd_bcf_reheader = "bcftools reheader --sample $file_new_header $old_vcfgz |bcftools view > $file_new_vcf" ;
                $log->andRun($cmd_bcf_reheader) ;

                my $file_orig = $old_vcfgz . ".orig" ;
                my $cmd_mvorig = "mv $old_vcfgz $file_orig" ;
                $log->andRun($cmd_mvorig) ;
                bgzip($file_new_vcf , "$file_new_vcf.gz") ;
                my $cmd_mvnewone = "mv $file_new_vcf.gz $old_vcfgz " ;
                $log->andRun($cmd_mvnewone) ;
            }
        }

    } else {
        $list_samples->{$sample} ++ ;
    }
}
close $SRC_ped ;
close $TGT_samplesex ;

my $file_ped_orig = $file_ped if($flag_multicallers) ;
my $file_sample_sex_orig = $file_sample_sex  if ($flag_multicallers) ;

if ($flag_multicallers) {
    close $TGT_pedmcaller ;
    close $TGT_samplesex_mcaller ;

    $file_ped = $file_ped_mcaller ;
    $file_sample_sex = $file_sample_sex_mcaller ;
}
$log->write ("Output sex infomation of Sample to $file_sample_sex") ;
$log->write ("Load ped file finish") ;

# InitAnnot 02: Input files preprocessing
#=======================================

my $log_stderr_vcf_concat = "$dir_log/stderr_vcfconcat.log" ;
my $log_stderr_vcf_sort = "$dir_log/stderr_vcfsort.log" ;
my $file_merged_vcfgz = "$file_prefix\_merged.vcf.gz" ;

my $files_to_process = [] ;
if (exists $sample2file->{all} ) {
    $log->write("found sample \"all\", skip concat and merge steps. ") ;

    my $file_allsample = shift @{$sample2file->{all}} ;
    if ($file_allsample =~ /\.gz$/) {
	$file_merged_vcfgz = $file_allsample ;

    } else {
	my $file_vcfgz_allsample = "$file_allsample.gz" ;
	my $cmd_bgzip_allsample = "bgzip $file_allsample -@ $threads_active -c > $file_vcfgz_allsample" ;
	$log->andRun($cmd_bgzip_allsample) ;
    }


} else {
    # InitAnnot 02-2: Concat vcf files from same sample 
    #=======================================
	      
    $log->write("Check and concat vcf files from same sample start.") ;

    foreach my $sample (keys %$sample2file) {
	my $filenum_in_one_sample = scalar @{$sample2file->{$sample}} ;
	if ($filenum_in_one_sample > 1) {
	    my $all_files = join (" ", @{$sample2file->{$sample}}) ;
	    my $file_concated_vcfgz = "$sample\_concated.vcf.gz" ;

	    my $cmd_vcf_concat = "vcf-concat $all_files 2>$log_stderr_vcf_concat | vcf-sort -p $threads_active 2>$log_stderr_vcf_sort | bgzip -@ $threads_active -c > $file_concated_vcfgz" ;
	    $log->andRun($cmd_vcf_concat) ;

	    push @$files_to_process , $file_concated_vcfgz ;

	} else {
	    push @$files_to_process , $sample2file->{$sample}->[0] ;
	}
    }

    foreach my $file_to_process (@$files_to_process) {
	tabix_vcf($file_to_process) ;
    }
    $log->write("Check and concat vcf files from same sample finish.") ;


    # InitAnnot 02-3: Merge vcf files from different sample
    #=======================================
    $log->write ("Merge vcf files from different samples start.") ;
    my $log_stderr_vcf_merge = "$dir_log/stderr_vcfmerge.log" ;

    my $cmd_merge = "vcf-merge -R '0/0' " . join (" ", @$files_to_process) . " 2> $log_stderr_vcf_merge |bgzip -@ $threads_active -c > $file_merged_vcfgz" ;
#    my $cmd_merge = "bcftools merge " . join (" ", @$files_to_process) . " 2> $log_stderr_vcf_merge |bgzip -@ $threads_active -c > $file_merged_vcfgz" ;
    $log->andRun($cmd_merge) ;
    $log->write ("Merge vcf files from different samples finish.") ;

    # InitAnnot 03: fix vcf ploidy
    #=======================================
    my $file_vcf_cleanXY = $file_prefix . "_cleanXY.vcf" ;
    my $file_vcfgz_cleanXY = $file_prefix . "_cleanXY.vcf.gz" ;
    my $log_stderr_fixploidy = "$dir_log/stderr_fixploidy.log" ;
    $log->write ("Fixing ploidy of sex chromosome start") ;

    my $cmd_fix_ploidy = "zcat $file_merged_vcfgz |vcf-fix-ploidy --samples $file_sample_sex >$file_vcf_cleanXY 2> $log_stderr_fixploidy" ;
    $log->andRun($cmd_fix_ploidy) ;

    my $cmd_bgzip_vcf_cleanXY = "bgzip -@ $threads_active -c $file_vcf_cleanXY > $file_vcfgz_cleanXY 2>> $log_stderr_fixploidy" ;
    $log->andRun("$cmd_bgzip_vcf_cleanXY") ;

    tabix_vcf($file_vcfgz_cleanXY) ;

    $log->write ("Fixing ploidy of sex chromosome finish") ;

	   
    # InitAnnot 05: VEP 
    #=======================================
    my $file_vcf_vep = $file_prefix . "_vep.vcf" ;
    my $file_vcfgz_vep = $file_prefix . "_vep.vcf.gz" ;
    my $script_runvep = "$dir_script/run_vep.pl" ;
    my $file_log_vep = "$dir_log/running_vep.$jobid.log" ;

    $log->write ("VEP start") ;
    my $cmd_run_vep = "$script_runvep -i $file_vcfgz_cleanXY -o $file_vcf_vep -t $jobid -n $threads_active -l $dir_log " ;
    $log->andRun($cmd_run_vep) ;
    $log->loadlog($file_log_vep) ;

    my $cmd_bgzip_vep = "bgzip -c -@ $threads_active $file_vcf_vep > $file_vcfgz_vep" ;
    $log->andRun($cmd_bgzip_vep) ;

    tabix_vcf($file_vcfgz_vep) ;
    $log->write ("VEP finish") ;

    # InitAnnot 06: Gemini load
    #=======================================
    $log->write ("Gemini load vcf to db start") ;

    print $CONFIG "DB\t$file_geminidb\n" ;
    print $CONFIG "VCF_ana\t$file_vcfgz_foranalysis\n" ;

    my $cmd_geminiload = "gemini load " ;
    $cmd_geminiload .= " --cores $threads_enable " ;
    $cmd_geminiload .= " -t VEP " ;
    $cmd_geminiload .= " -v $file_vcfgz_vep " ;
    $cmd_geminiload .= " -p $file_ped " ;
    $cmd_geminiload .= " $file_geminidb 2> $dir_log/stderr_geminiload_$jobid.log" ;


    $log->andRun($cmd_geminiload) ;
    $log->write("Gemini load vcf to db finish") ;

    my $cmd_cp_foranalysis = "cp $file_vcfgz_cleanXY $file_vcfgz_foranalysis" ;
    $log->andRun($cmd_cp_foranalysis) ;


    # snpEff
    #=======================================
    my $file_vcfgz_snpeff = $file_prefix . "_snpeff.vcf.gz" ;
    my $script_run_snpeff = "$dir_script/run_snpeff.pl" ;
    my $file_log_snpeff = "$dir_log/running_snpeff.$jobid.log" ;

    $log->write("Run snpEff start") ;
    my $cmd_run_snpeff = "$script_run_snpeff -i $file_vcfgz_cleanXY -o $file_vcfgz_snpeff -j $jobid -n $threads_active -d $file_geminidb" ;
    $log->andRun($cmd_run_snpeff) ;
    $log->write("Run snpEff finish") ;

    # ANNOVAR
    #=======================================
    $log->write("Running ANNOVAR start") ;
    my $file_vcfgz_annovar = "$file_prefix\_annovar.vcf.gz" ;
    my $file_log_running_annovar = "$dir_log/running_annovar.log" ;
    my $cmd_runannovar = "$dir_script/run_annovar.pl" ;
#    $cmd_runannovar .= " -i $file_vcfgz_cleanXY " ;
    $cmd_runannovar .= " -i $file_vcfgz_snpeff " ;
    $cmd_runannovar .= " -o $file_vcfgz_annovar " ;
    $cmd_runannovar .= " -j $jobid " ;

    $log->write($cmd_runannovar) ;
    `$cmd_runannovar` ;

    tabix_vcf($file_vcfgz_annovar) ;
    $log->write("Running ANNOVAR finish") ;

    # Annotate ANNOVAR
    #=======================================
    $log->write("Extract info fields from ANNOVAR vcf start") ;
    my $annovar_header_parser = Vcf_parser->new (file => $file_vcfgz_annovar) ;

    my $list_extract = [] ;
    my $list_type = [] ;
    my $list_column = [] ;
    my $list_operation = [] ;

    my $gannot_list = retrieve $file_hashref_vannot_tmp ;

    foreach my $info_tag (sort keys %{$annovar_header_parser->{header}->{INFO}} ) {
        next if ($info_tag eq 'lines') ;  # dirty hack for Vcf_parser bug

        my $desc = $annovar_header_parser->{header}->{INFO}->{$info_tag}->{Description} ;

        if ($desc =~ /provided by ANNOVAR/) {
	    my $col_tag = $info_tag ;
	    $col_tag =~ s/\-/\_/g ;
	    $col_tag =~ s/\+/x/g ;

	    next if ($info_tag eq 'esp6500siv2_all') ; # dirty hack for duplicate INFO tag, there is ESP6500siv2_ALL already

	    push @{$gannot_list->{extract_list}}   , $info_tag ;
	    push @{$gannot_list->{type_list}}      , 'text' ;
	    push @{$gannot_list->{coladd_list}}    , "anv_$col_tag" ;
	    push @{$gannot_list->{operation_list}} , 'list' ;

#            push @$list_extract , $info_tag ;
#            push @$list_type , 'text' ;
#            push @$list_column , "anv_$col_tag" ;
#            push @$list_operation , 'list' ;
        }
    }

#    my $cmd_gannot_annovar = "gemini annotate -f $file_vcfgz_annovar -a extract " ;
#    $cmd_gannot_annovar .= " -e " . join("," , @$list_extract) ;
#    $cmd_gannot_annovar .= " -t " . join("," , @$list_type) ;
#    $cmd_gannot_annovar .= " -c " . join("," , map {"'" . $_ . "'"} @$list_column) ;
#    $cmd_gannot_annovar .= " -o " . join("," , @$list_operation) ;
#    $cmd_gannot_annovar .= " $file_geminidb 2> $dir_log/stderr.gannot_annovar.$jobid.log" ;
#    $log->andRun($cmd_gannot_annovar) ;

    $log->write("Extract info fields from ANNOVAR finish") ;


    # Varapp Support (Gemini annotate with AF,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR columns)
    #=======================================
    $log->write("Extract info field from vcf for Varapp (AF,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR) start") ;
    my $cmd_gemini_annotate = "gemini annotate -f $file_vcfgz_vep -a extract " ;
    $cmd_gemini_annotate .= " -e AF,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR " ; # extract from INFO, for Varapp compatiable
    $cmd_gemini_annotate .= " -t float,float,float,float,float,float " ;
    $cmd_gemini_annotate .= " -c AF,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR " ;
    $cmd_gemini_annotate .= " -o mean,mean,mean,mean,mean,mean " ;
    $cmd_gemini_annotate .= "  $file_geminidb 2> $dir_log/stderr.gannot_forvarapp.$jobid.log" ;

#    $log->andRun($cmd_gemini_annotate) ;

    my $varapp_required_list = [qw/AF BaseQRankSum FS MQRankSum ReadPosRankSum SOR/] ;
    push @{$gannot_list->{extract_list}}   , @$varapp_required_list ;
    push @{$gannot_list->{type_list}}      , map {'float'} @$varapp_required_list ;
    push @{$gannot_list->{coladd_list}}    , @$varapp_required_list ;
    push @{$gannot_list->{operation_list}} , map {'mean'} @$varapp_required_list ;

    $log->write("Extract info field from vcf for Varapp (AF,BaseQRankSum,FS,MQRankSum,ReadPosRankSum,SOR) finish") ;



    # Extract pathway info
    #=======================================
    $log->write("gemini pathway start") ;

    my $file_vcfgz_pathway = "$file_prefix\_pathway.vcf.gz" ;

    my $cmd_running_pathways = "$dir_script/pathway_parse.pl " ;
    $cmd_running_pathways   .= " -j $jobid " ;
    $cmd_running_pathways   .= " -d $file_geminidb " ;
#    $cmd_running_pathways   .= " -v $file_vcfgz_cleanXY " ;
    $cmd_running_pathways   .= " -v $file_vcfgz_annovar " ;
    $cmd_running_pathways   .= " -l $dir_log " ;
    $cmd_running_pathways   .= " -o $file_vcfgz_pathway " ;
    $log->andRun($cmd_running_pathways) ;
    $log->write("gemini pathway finish") ;

    push @{$gannot_list->{extract_list}}   , 'pathway' ;
    push @{$gannot_list->{type_list}}      , 'text' ;
    push @{$gannot_list->{coladd_list}}    , 'pathway' ;
    push @{$gannot_list->{operation_list}} , 'list' ;

    # Extract GO infomation
    #=======================================
    $log->write("GO parse start") ;

    my $file_vcfgz_GO = "$file_prefix\_go.vcf.gz " ;

    my $cmd_parse_GO = "$dir_script/GO_parse.pl " ;
    $cmd_parse_GO   .= " -j $jobid " ;
    $cmd_parse_GO   .= " -d $file_geminidb " ;
#    $cmd_parse_GO   .= " -v $file_vcfgz_cleanXY " ;
    $cmd_parse_GO   .= " -v $file_vcfgz_pathway " ;
    $cmd_parse_GO   .= " -l $dir_log " ;
    $cmd_parse_GO   .= " -o $file_vcfgz_GO " ;

    $log->andRun($cmd_parse_GO) ;
    $log->write("GO parse finish") ;

    my $file_log_goparse = "$dir_log/log_goparse_$jobid.log"  ;
    $log->loadlog($file_log_goparse) ;
   
    my $GO_cols = [qw/GO_Biological_Process GO_Cellular_Component GO_Molecular_Function/] ;
    push @{$gannot_list->{extract_list}}   ,  @$GO_cols ;
    push @{$gannot_list->{type_list}}      , map {'text'} @$GO_cols ;
    push @{$gannot_list->{coladd_list}}    , @$GO_cols  ;
    push @{$gannot_list->{operation_list}} , map {'list'} @$GO_cols ;


    store $gannot_list , $file_hashref_vannot_tmp ;

    my $cmd_gannot = "$dir_script/gannot.pl -j $jobid -d $file_geminidb -g $file_hashref_vannot_tmp -v $file_vcfgz_GO " ;
    $log->andRun($cmd_gannot) ;
}
               
close $CONFIG ;

my $cmd_cp_config = "cp $file_config $file_config_db" ;
$log->andRun($cmd_cp_config) ;

$log->write("VSannot finish") ;


sub tabix_vcf {
    my $file = shift ;

    if ($file !~ /gz$/) {
        my $file_orig = $file ;
        $file = "$file_orig.gz" ;
        bgzip($file_orig,$file) ;
    }

    my $cmd_tabix_vcf = "tabix -p vcf -f $file " ;
    $log->andRun($cmd_tabix_vcf) ;
}

sub bgzip {
    my $file_src = shift ;
    my $file_tgt = shift ;

    my $cmd_bgzip = "bgzip -@ $threads_active -c $file_src > $file_tgt" ;
    $log->andRun($cmd_bgzip) ;
}


