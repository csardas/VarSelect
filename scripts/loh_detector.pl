#!/usr/bin/perl

# detect LOH in paired samples

use strict ;
use Getopt::Std ;
use Data::Dumper ;
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';
use VarSelect ;
use VarSelect::Log ;
use VarSelect::Vcf_parser ;


# defaule global setting
my $dir_script = dirname(abs_path $0) ;
my $dir_vsHOME = dirname($dir_script) ;
my $dir_db = "$dir_vsHOME/db" ;
my $dir_lib = "$dir_vsHOME/lib" ;
my $dir_workflows = "$dir_vsHOME/workflows" ;

# Options handle
my $opts = {} ;
getopts("j:i:o:p:l:d:n:h", $opts) ;

die ( usage() ) if $opts->{h} ;
my $jobid	= $opts->{j} || getts() ;
my $dir_log	= $opts->{l} ;
my $file_ped	= $opts->{p} || die "Error: -p PED file is required!\n\n" . usage() ;
my $file_vcf	= $opts->{i} || die "Error: -v VCF file is required!\n\n" . usage() ;
my $file_db	= $opts->{d} ;
my $threads_enable = $opts->{n} || 16 ;
my ($input_name,$input_path,$input_ext) = fileparse($file_vcf, qw/.vcf .vcf.gz/) ;

my $file_output_vcf = $opts->{o} || $input_path . $input_name . "\_loh.vcf.gz" ;

# jobs specific path
my $dir_results_output = "./VarSelectAnalysisResult_$jobid" ;
my $dir_log = "$dir_results_output/log_$jobid" ;
my $dir_work = "$dir_results_output/work_$jobid" ;

die "Current directory is not writable! Please check your permission!\n" unless (-w "./" ) ;
make_path($dir_results_output , { chmod => 0755 ,} ) unless (-e $dir_results_output) ;
make_path($dir_log , $dir_work ,{ chmod => 0755 ,} ) unless (-e $dir_work && -e $dir_log) ;

my $file_output_isloh = "$dir_results_output/$input_name\_$jobid\_loh.txt" ;

# start loging
my $file_log = "$dir_log/log_lohdetect_$jobid.log" ;
my $log = VarSelect::Log->new(file => $file_log) ;

$log->write("LOH detecting start") ;

# Step 1. Load PED sample info
open (my $SRC_ped , "$file_ped") || die "Can't open $file_ped!\n$!\n" . usage() ;

my $affect_sample_list = [] ;
while (<$SRC_ped>) {
    chomp ;
    next if /^#/ ;
    my ($family,$sample,$father,$mother,$sex,$phenotype) = split /\t/ ;

    push @$affect_sample_list , $sample if ($phenotype == 2) ;
}
close $SRC_ped ;
my $affect_sample_hash = {map {$_ => 1} @$affect_sample_list} ;

# Step 2. Load VCF and output detect result
my $vcf_parser = Vcf_parser->new (file => $file_vcf) ;
open (my $TGT_rst , ">$file_output_isloh") ;
print $TGT_rst  join("\t",(qw/#CHROM POS REF ALT is_LOH/, @{$vcf_parser->{samples}}) ) . "\n" ;

while (my $var = $vcf_parser->next_var) {
    my $is_LOH = 0 ;

    # Step 3. check GT
    my $flag_have_homo_affected_sample = 0 ;
    my $homo_affected_gt = [] ;
    my $affected_gt = [] ;
    my $unaffected_gt = [] ;

    foreach my $sample (@{$var->{samples}} ) {
	my $genotype = $var->{sample_val}->{$sample}->{GT} ;

	$genotype =~ s/\./0/g ; # treat . as 0
	if (sample_is_affected($sample)) {
	    if (genotype_is_homo($genotype) and !genotype_is_refref($genotype) ) {	# GT1 get
		$flag_have_homo_affected_sample = 1 ;
		push @$homo_affected_gt , $genotype ;
	    }
	    push @$affected_gt , $genotype ;
	} else {
	    push @$unaffected_gt , $genotype ;
	}
    }

    my $unaff_gt_hash = array2hash($unaffected_gt);

    my $flag_have_same_gt_in_homo_aff_and_unaff ;

    if ($flag_have_homo_affected_sample) {

	foreach my $hm_af_gt (@$homo_affected_gt) {
	    $flag_have_same_gt_in_homo_aff_and_unaff = 1 if (exists $unaff_gt_hash->{$hm_af_gt} ) ;
	}

	unless ($flag_have_same_gt_in_homo_aff_and_unaff) {
	    foreach my $ua_gt (@$unaffected_gt) {
		unless (genotype_is_homo($ua_gt)) { # GT2 get

		    my ($alleltype1, $alleltype2) = split /[\|\/]/ , $ua_gt ;
		    foreach my $h_a_gt (@$homo_affected_gt) {
			my ($h_allele1,$h_allele2) = split /[\|\/]/ , $h_a_gt ;

			if ($h_allele1 eq $alleltype1 || $h_allele1 eq $alleltype2) { # GT1 X GT2
			    $is_LOH = 1 ;

			}
		    }
	    
		}
	    }
	}
    }
    # Step 4. output is_LOH 
    print $TGT_rst join("\t" , map {$var->{$_} } qw/CHROM POS REF ALT/ ) ;
    print $TGT_rst "\t$is_LOH" ;

    foreach my $sample (@{$var->{samples}}) {
	print $TGT_rst "\t$sample:" . $var->{sample_val}->{$sample}->{GT} ;
    }
    print $TGT_rst "\n" ;

}

# Step 5. vcf-annotate
my $cmd_bgzip_output = "bgzip -@ $threads_enable $file_output_isloh" ;
$log->andRun($cmd_bgzip_output) ;

my $cmd_tabix_output = "tabix -f -s 1 -b 2 -e 2 -c # $file_output_isloh.gz " ;
$log->andRun($cmd_tabix_output) ;

#my $cmd_create_new_vcf = "cp $file_vcf $file_output_vcf" ;
#$log->andRun($cmd_create_new_vcf) ;

my $cmd_vcfannotate = "zcat $file_vcf | vcf-annotate -a $file_output_isloh.gz " ;
$cmd_vcfannotate .= " -c " . join("," , qw"CHROM POS REF ALT INFO/is_LOH -") ;
$cmd_vcfannotate .= " -d key=INFO,ID=is_LOH,Number=1,Type=Integer,Description='is_LOH annotation provided by Varnotate' " ;
$cmd_vcfannotate .= " 2> $dir_log/stderr_vcfannotate_$jobid.log" ;
$cmd_vcfannotate .= " | bgzip -c > $file_output_vcf " ;
$log->andRun($cmd_vcfannotate) ;

my $cmd_tabix_new_vcf = "tabix -p vcf $file_output_vcf" ;
$log->andRun($cmd_tabix_new_vcf) ;

# Step 6. gemini annotate
my $cmd_gannotate_LOH = "gemini annotate -f $file_output_vcf " ;
$cmd_gannotate_LOH .= " -a extract " ;
$cmd_gannotate_LOH .= " -e is_LOH " ;
$cmd_gannotate_LOH .= " -t integer " ;
$cmd_gannotate_LOH .= " -c is_LOH_$jobid " ;
$cmd_gannotate_LOH .= " -o first " ;
$cmd_gannotate_LOH .= " $file_db 2> $dir_log/stderr.geminiannotate_loh.$jobid.log" ;

$log->andRun($cmd_gannotate_LOH) ;

$log->write("LOH detecting finish") ;

sub array2hash {
    my $input = shift ;
    return { map {$_=>1} @$input } ;
}

sub genotype_is_homo {
    my $gt = shift ;
    my @alleletype = split /[\/\|]/ , $gt ;

    return ($alleletype[0] eq $alleletype[1])?1:0 ;
}

sub genotype_is_refref {
    my $gt = shift ;
    my @alleletype = split /[\/\|]/ , $gt ;

    return ($alleletype[0] eq '0' && genotype_is_homo($gt) )?1:0 ;
}

sub sample_is_affected {
    my $sample = shift ;
    return (exists $affect_sample_hash->{$sample})?1:0 ;
}

sub usage {
return <<EOUsage
$0 version 
Usage:
    $0 -i file -o file -p file

    -i  VCF file input                      (Required)
    -o  VCF file output                     (Required)
    -p  PED file                            (Required)

    -d  gemini db
    -n  number of threads

    -h  this page

EOUsage
}

