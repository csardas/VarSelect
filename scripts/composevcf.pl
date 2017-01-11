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

# options handle
my $opts = {} ;
getopts("v:h", $opts) ;
my $jobid = $opts->{j} ;
my $file_input = $opts->{v} ;
my $vcf_parser = Vcf_parser->new (file => $file_input) ;

my $prevar ;
my $count = {} ;
while (my $var = $vcf_parser->next_var) {
    $count->{all} ++ ;

    if ($prevar && $var->{POS} == $prevar->{POS}) {
	print "prevar:\t$prevar->{CHROM}\t$prevar->{POS}\t$prevar->{REF}\t$prevar->{ALT}\n" ;
	print "var:\t$var->{CHROM}\t$var->{POS}\t$var->{REF}\t$var->{ALT}\n" ;

	my $len_p_ref = length($prevar->{REF}) ;
	my $len_v_ref = length($var->{REF}) ;
	my $len_p_alt = length($prevar->{ALT}) ;
	my $len_v_alt = length($var->{ALT}) ;

	my $sample = $var->{samples}->[0] ;

	$count->{decompose} ++ ;

	if ($len_v_ref == 1 && $len_p_ref == 1) {
	    if ($prevar->{REF} eq $var->{REF}) {
		$count->{simple} ++ ;
		print "simple:\t$var->{REF}\t$prevar->{ALT},$var->{ALT}\n\n" ;
	    } else {
		$count->{complex} ++ ;
	    }
	} elsif ($len_v_ref > $len_p_ref) {

	    if ($len_p_alt == $len_v_alt && $len_v_alt == 1) {

		if ($var->{ALT} eq $prevar->{ALT}) {
		    $count->{del_v} ++ ;
		    my $alt2 = $var->{ALT} ;
		    my $alt1patt = substr ($prevar->{REF},1) ;

		    my $alt1 = $var->{REF} ;
		    $alt1 =~ s/$alt1patt// ;

		    print "delV:$var->{REF}\t$alt1,$alt2\n(alt1patt = $alt1patt)\n\n" ;
		} else {
		    print "diffalt!\n\n" ;
		    $count->{complex} ++ ;
		}
	    } else {
		$count->{complex} ++ ;
		print "Complex!\n\n" ;
	    }


	    print "$var->{REF}\n\n" ;
	} elsif ($len_v_ref < $len_p_ref) {

	    if ($len_p_alt == $len_v_alt && $len_v_alt == 1) {


		if ($var->{ALT} eq $prevar->{ALT}) {
		    $count->{del_p} ++ ;

		    my $alt1 = $prevar->{ALT} ;
		    my $alt2patt = substr($var->{REF},1) ;
		    my $alt2 = $prevar->{REF} ;
		    $alt2 =~ s/$alt2patt// ;

		    print "delpa=:$prevar->{REF}\t$alt1,$alt2\n\n" ;
		} else {
		    $count->{complex} ++ ;
		}
	    } else {
		$count->{complex} ++ ;
		print "Complex!\n\n" ;
	    }

	    print "$prevar->{REF}\n\n" ;
	} else {
	    $count->{other} ++ ;
	    print "other\n\n" ;
	}
	print "INFO: $var->{INFO}->{all}\n" ;
	print "INFO: $prevar->{INFO}->{all}\n" ;

	print "FORMAT: $var->{FORMAT}\n" ;
	print "FORMAT: $prevar->{FORMAT}\n" ;

	my $sample = $var->{samples}->[0] ;
	print "Sample: $var->{$sample}\n" ;
	print "Sample: $prevar->{$sample}\n" ;

	print "Sample: " . combine_vcf_sample($prevar,$var) . "\n" ;
    }

    

    $prevar = $var ;
}

foreach my $idx (sort keys %$count) {
    print "$idx\t$count->{$idx}\n" ;
}

sub combine_vcf_sample {
    my $var1 = shift ;
    my $var2 = shift ;

    my $sampleid = $var1->{samples}->[0] ;
    my $sample1 = $var1->{$sampleid} ;
    my $sample2 = $var2->{$sampleid} ;

    my @array1 = split(/:/,$sample1) ;
    shift @array1 ;
    my $sample1_nogt = join(":" , @array1) ;

    my @array2 = split(/:/,$sample2) ;
    shift @array2 ;
    my $sample2_nogt = join(":" , @array2) ;

    if ($sample1_nogt eq $sample2_nogt) {
	return "1/2:$sample1_nogt" ;
    } else {
	return "1/2" ;
    }

}

sub is_heterozygous {
}

