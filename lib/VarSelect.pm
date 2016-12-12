package VarSelect ;
use Exporter () ;
@ISA = qw(Exporter);
@EXPORT= qw(getts tabix_vcf bgzip $Setting);
use vars qw($Setting) ;

$Setting = {
     dir_vep_cache  => "/peterhouse/tools/bcbio/genomes/Hsapiens/GRCh37/vep" ,
     dir_vep	    => "/peterhouse/tools/bcbio/anaconda/bin" ,       
     vep	    => "variant_effect_predictor.pl" ,
     file_ref	    => "/peterhouse/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa" ,
     file_dbNSFP    => "/peterhouse/lchen/csardas/varselect/db/dbNSFP_3.2a.txt.gz" ,
     jar_snpeff	    => "/peterhouse/tools/bcbio/anaconda/share/snpeff-4.3-3/snpEff.jar" ,
     dir_snpeff	    => "/peterhouse/tools/bcbio/genomes/Hsapiens/GRCh37/snpeff" ,
     genome_snpeff  => "GRCh37.75" ,
     dir_annovar    => "/peterhouse/tools/annovar/" ,
     genome_annovar => "hg19" ,
} ;

sub getts {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    return sprintf ( "%04d%02d%02d%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
}


sub tabix_vcf {
    my $file = shift ;
    my $vs_log = shift || 0 ;

    if ($file !~ /gz$/) {
        my $file_orig = $file ;
        $file = "$file_orig.gz" ;
        bgzip($file_orig,$file,$vs_log) ;
    }

    my $cmd_tabix_vcf = "tabix -p vcf -f $file " ;

    if ($vs_log) {
	$vs_log->andRun($cmd_tabix_vcf) ;
    } else {
	`$cmd_tabix_vcf` ;
    }
}

sub bgzip {
    my $file_src = shift ;
    my $file_tgt = shift ;
    my $vs_log = shift || 0 ;
    my $cores = shift || 16 ;

    my $cmd_bgzip = "bgzip -@ $cores -c $file_src > $file_tgt" ;

    if ($vs_log) {
	$vs_log->andRun($cmd_bgzip) ;
    } else {
	`$cmd_bgzip`
    }
}

1;
