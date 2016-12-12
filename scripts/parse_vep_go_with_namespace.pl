#!/usr/bin/perl

use strict ;
use Getopt::Std ;
use Storable ;
use File::Basename ;
use Cwd  qw(abs_path);
use lib dirname(dirname( abs_path($0) )) . '/lib';

use VarSelect ;
use VarSelect::Vcf_parser ;

my $dir_script = dirname(abs_path $0) ;
my $dir_vs_home = dirname($dir_script) ;

print "dir_script:$dir_script\n" ;
print "vshome: $dir_VarSelect_home\n" ;

my $opts = {} ;
getopts("v:o:h") ;
die usage() if($opts->{h}) ;
my $file_vcfgz_input = $opts->{v} || die "\n\n-v VCF_file is required!\n\n" . usage() ;


sub usage {
return <<EOUsage
$0 

VarSelect Component

EOUsage
}

