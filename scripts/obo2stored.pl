#!/usr/bin/perl 
use strict ;
use Storable ;
use YAML ;

my $obofile = "./gene_ontology.1_2.obo" ;
my $dir_dot = "dots" ;
my $dir_png = "graphics" ;
my $file_go_stored = "./goterm_1_2.hashref" ;

my $data ;
my $n = 0 ;
my $GO = {} ;
my $isa = {} ;



# read obo file
open (SRC,$obofile) ;
while (<SRC>) {
    $data .= $_ ;

    if (($_ eq "[Term]\n") || eof(SRC) ) {

	my $term = parser_Term($data) ;
	$GO->{$term->{id}} = $term ;

	$data = '' ;
    }
}
close SRC ;

print "reading done...." ;

store $GO , "$file_go_stored" ;



sub parser_Term {
    my $data = shift ;
    my $output = {} ;

#[Term]
#id: GO:0000001
#name: mitochondrion inheritance
#namespace: biological_process
#def: "The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824, PMID:11389764]
#synonym: "mitochondrial inheritance" EXACT []
#is_a: GO:0048308 ! organelle inheritance
#is_a: GO:0048311 ! mitochondrion distribution
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
