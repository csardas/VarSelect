#!/usr/bin/perl

use strict ;
use Storable ;

my $mapping_ENS2GO = "../db/ENS2GO.txt" ;
my $mapping_go2ns = "../db/go2ns.hashref" ;
my $file_hashref_enst2go = "../db/ENST2GO.hashref" ;


my $ns2 = {
    biological_process => 'BP' ,
    cellular_component => 'CC' ,
    molecular_function => 'MF' ,
} ;

my $GO2NameSpace = retrieve "$mapping_go2ns" ;

open (my $ENS2GO , "$mapping_ENS2GO") ;

my $header = <$ENS2GO> ;
chomp $header ;


# my $header_col = [split /\t/, $header] ;
# Ensembl Gene ID        Ensembl Transcript ID   Ensembl Protein ID      GO Term Accession       GO Term Name

my $data = {} ;
my $nslist = {} ;
while (my $line = <$ENS2GO>) {
    chomp $line ;
    my ($ensg,$enst,$ensp,$gotermid,$gotermname) = split /\t/,$line ;

    my $namespace = $GO2NameSpace->{$gotermid} ;
    $nslist->{$namespace} = 1 ; 
    $gotermname =~ s/ /\_/g ;

    if ($gotermid) {
	push @{$data->{$enst}->{$namespace}} , "$gotermid:$gotermname" ;    
    }
}

my $output = {} ;
foreach my $transcript (sort keys %$data) {
    print "$transcript" ; 
    foreach my $namespace (sort keys %$nslist ) {

	next unless $namespace ;

	if (exists $data->{$transcript}->{$namespace} ) {
	    my $go_text = join("&" , @{$data->{$transcript}->{$namespace}}) ;
	    print "\t$go_text" ;
	    $output->{$transcript}->{$ns2->{$namespace}} = $go_text ;
	} else {
	    $output->{$transcript}->{$ns2->{$namespace}} = '' ;
	    print "\t" ;
	}
    }
    print  "\n";
#map { join("," , @{$data->{$transcript}->{$_}} ) } sort keys %$nslist) . "\n" ;
}

store $output , $file_hashref_enst2go ;
