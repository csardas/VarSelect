#!/usr/bin/perl

our $DEFAULT_ref_GRCh37 = "/peterhouse/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa" ;
our $dir_vep_cache = "/peterhouse/tools/bcbio/genomes/Hsapiens/GRCh37/vep" ;

our $dir_vep = "/peterhouse/tools/bcbio/anaconda/bin" ;       # PATH setting for VEP on Sydney
our $exe_vep = "variant_effect_predictor.pl" ;
our $file_ref = $DEFAULT_ref_GRCh37 ;
our $file_dbNSFP = "/peterhouse/lchen/csardas/varselect/db/dbNSFP_3.2a.txt.gz" ;

1;
