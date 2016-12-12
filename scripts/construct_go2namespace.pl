#!/usr/bin/perl

use strict ;
use Data::Dumper ;
use Storable ;

my $url_go_basic = "http://geneontology.org/ontology/go-basic.obo" ;

my $file_obo_go_basic = "../db/go-basic.obo" ;
my $file_hashref_go2ns = "../db/go2ns.hashref" ;

my $go2ns = {} ;
my $data ;
open (my $SRC, $file_obo_go_basic) ;
while (<$SRC>) {
    $data .= $_ ;

    if (($_ eq "[Term]\n") || eof($SRC) ) {

        my $term = parser_Term($data) ;
	$go2ns->{$term->{id}} = $term->{namespace} ;

        $data = '' ;
    }
}
close $SRC ;

#print Dumper $go2ns ;

store $go2ns , $file_hashref_go2ns ;


sub parser_Term {
    my $data = shift ;
    my $output = {} ;

    my @lines = split(/\n/,$data) ;

    my $isa_array = [] ;
    for (@lines) {

        my ($idx,$val) = split(/\: /,$_,2) ;

        if ($idx eq "is_a") {
            my ($goid,$goname) = split(/ \! /,$val) ;
            push @$isa_array,$goid ;
        } else {
            $output->{$idx} = $val ;
        }

    }
    $output->{is_a} = $isa_array ;

    return $output ;
}

#[Term]
#id: GO:0000001
#name: mitochondrion inheritance
#namespace: biological_process
#def: "The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824, PMID:11389764]
#synonym: "mitochondrial inheritance" EXACT []
#is_a: GO:0048308 ! organelle inheritance
#is_a: GO:0048311 ! mitochondrion distribution

