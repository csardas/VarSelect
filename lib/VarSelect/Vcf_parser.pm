package Vcf_parser ;

use strict ;
use Carp ;

sub new {
    my $invocant = shift ;
    my $class	 = ref $invocant || $invocant ;
    my %initial_attritubes = @_ ;

    my $self = {
        file	=> undef ,
        fh	=> undef ,
        string	=> undef ,
        %initial_attritubes ,
    } ;

    # Attributes: 
    #	count_varnum_first

    my $all_varnum = 0 ;
    my $flag_count_varnum_first = $self->{count_varnum_first} || 0 ;

    my $SRC ;
    if ($self->{file} =~ /\.bz2$/) { 
        open ($SRC, "bzcat $self->{file}|") or confess "Error: Couldn't open file $self->{file} or can't use bzcat to decompress $self->{file} .\n" ;

	$all_varnum = `bzgrep -c -v '^#' $self->{file} ` if ($flag_count_varnum_first) ;

    } elsif($self->{file} =~ /\.gz$/) {
	open ($SRC, "zcat $self->{file}|") or confess "Error: Couldn't open file $self->{file} or can't use zcat to decompress $self->{file} .\n" ;

	$all_varnum = `zgrep -c -v '^#' $self->{file} ` if ($flag_count_varnum_first) ;

    } elsif(defined $self->{file}) {
        open ($SRC, "$self->{file}") or confess "Error: Couldn't open file $self->{file}\n" ;

	$all_varnum = `grep -c -v '^#' $self->{file}` if ($flag_count_varnum_first) ;

    } elsif (defined $self->{fh}) {
        $SRC = $self->{fh} ;
	$self->{count_varnum_first} = 0 ;

    } elsif (defined $self->{string}) {
        open $SRC , "<" , \{$self->{string}} ;
	$self->{count_varnum_first} = 0 ;

    } else {
        croak "option 'file' or 'fh' or 'string' should be defined\n" ;
    }

    chomp $all_varnum if($flag_count_varnum_first) ;

    my $header = {} ;
    my $headerline = '' ;
    my $headerline_num = 0 ;
    my $colname = [] ;
    my $colline = '' ;
    my $samples = [] ;

    my $fixed_fields = { map {$_ => 1} qw/CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/} ;

    while (my $linein = <$SRC>) {
	chomp $linein ;
        my $char1 = substr($linein,0,1) ;

        if ($char1 eq '#') {	 # start with # 
	    my $char2 = substr($linein,1,1) ;

	    if ($char2 eq '#') { # start with ##  => header line
		my $line = substr($linein,2) ;

		my ($h_key,$h_val) = split(/\=/, $line , 2 ) ;
		push @{$header->{$h_key}->{lines}} , $h_val ;

		$headerline .= "$linein\n" ;
		$headerline_num ++ ;

	    } else { # start with single #  => column name line
		my $line = substr($linein,1) ;
		@$colname = split (/\t/, $line) ;
		$colline = $linein ;
		$headerline .= "$linein\n" ;
		$headerline_num ++ ;
		foreach my $col (@$colname) {
		    push @$samples , $col unless (exists $fixed_fields->{$col}) ;
		}

		last ;
	    }
	}
    }

    foreach my $hkey (keys %$header) {
	foreach my $headerline (@{$header->{$hkey}->{lines}}) {
	    my $hash = _parse_header_line($hkey,$headerline) ;

	    if (exists $hash->{ID}) {
		$header->{$hkey}->{$hash->{ID}} = $hash ;
	    } else {
		$header->{$hkey}->{value} = $hash->{value} ;
	    }
	}
    }

    $self->{header} = $header ;
    $self->{headerline} = $headerline ;
    $self->{colname} = $colname ;
    $self->{colnum} = $#$colname + 1 ;
    $self->{colline} = $colline ;
    $self->{file_handle} = $SRC ;
    $self->{samples} = $samples ;
    $self->{all_varnum} = $all_varnum if ($flag_count_varnum_first) ;

    bless ($self,$class) ;
    return $self ;

}

sub next_var {
    my $self = shift ;
    my $linedata = {} ;
    local $/ = "\n" ;
    my $samples = $self->{samples} ;


    my $fh = $self->{file_handle} ;
    return unless my $linein = <$fh> ;
    chomp $linein ;
    
    my @data = split(/\t/,$linein) ;
    
    for (my $i = 0 ; $i < $self->{colnum} ; $i++ ) {
        $linedata->{ $self->{colname}->[$i] } = $data[$i] ;
    }
    my @formats = split(/\:/, $linedata->{FORMAT}) ;
    

    $linedata->{line} = $linein ;

    foreach my $sample ( grep {$linedata->{$_} ne '.' }  @$samples ) {
	my @data = split(/\:/,$linedata->{$sample} ) ;
	%{$linedata->{sample_val}->{$sample}} = map {$formats[$_] => $data[$_] } 0 .. $#formats ;
    }

    my $INFO = {} ;
#    %$INFO = map { split(/\=/,$_) } split(/\;/, $linedata->{INFO}) ;
    %$INFO = map {(/\=/)?split (/\=/,$_,2):($_,'True')} split(/\;/, $linedata->{INFO}) ;
    $INFO->{all} = $linedata->{INFO} ;

    $linedata->{INFO} = $INFO ;
    $linedata->{samples} = $samples ;

    bless $linedata , "Vcf" ;
    return $linedata ;
}

sub header {
    my $self = shift ;
    return $self->{headerline} ;
}

sub _parse_header_line {
    # this function is inferred from "sub parse_header_line {}" of Vcf.pm from vcftools
    my $key = shift ;
    my $value = shift ;
#    if ( !($line=~/^([^=]+)=/) ) { return { key=>$line, value=>'' }; }
#    my $key   = $1;
#    my $value = $';

    my $desc;
    if ( $value=~/\s*\"([^\"]+)\"\s*>$/ ) {
	$desc  = $1 ; 
	$value = $` ; 
    }

    return { key=>$key, value=>$value } unless ($desc) ;

    if ( $key eq 'INFO' or $key eq 'FORMAT' ) {
        my ($id,$number,$type,@rest) = map {(split(/\=/,$_))[1]} split(/,\s*/,$value);
#        if ( !$type or scalar @rest ) { $self->throw("Could not parse the header line: $line\n"); }
        return { key=>$key, ID=>$id, Number=>$number, Type=>$type, Description=>$desc };
    }

    if ( $key eq 'FILTER' or $key eq 'ALT') {
        my ($id,@rest) = map {(split(/\=/,$_))[1]} split(/,\s*/,$value);
#        if ( !$id or scalar @rest ) { $self->throw("Could not parse the header line: $line\n"); }
        return { key=>$key, ID=>$id, Description=>$desc };
    }
#    $self->throw("Could not parse the header line: $line\n");
#
#   Other header record
    my ($id,@rest) = map {(split(/\=/,$_))[1]} split(/,\s*/,$value);
    return { key=>$key, ID=>$id, Description=>$desc };
}

sub DESTROY {
    my $self = shift ;
    close $self->{file_handle} ;
}

1;

#self
#   ->{header}
#	->{INFO}->{lines}->[0] = "<ID=MDV,Number=1,Type=Integer,Description="Maximum number of high-quality nonRef reads in samples">"
#		         ->[1] = "<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias (v2) for filtering splice-site artefacts in RNA-seq data. Note: this version may be broken.">"
#			 ->[2] = <ID=SF,Number=.,Type=String,Description="Source File (index to sourceFiles, f when filtered)">
#			.....
#			.....
#
#	->{FORMAT}->{lines}->[0] = <ID=GT,Number=1,Type=String,Description="Genotype">
#			    ->[1] = <ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#			    ->[2] = <ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
#			    .......
#			    .......
#
#   ->{samples} = [W164, W173, W180-3, W182]
#
#   ->{INFO}->{col} = "AC1=2;AC=2;AF1=1;AN=2;DP4=0,0,2,0;DP=2;FQ=-33;MQ=20;SF=26;VDB=7.520000e-02"
#	    ->{AC1} = 2
#	    ->{AC}  = 2
#	    ->{AF1} = 1
#	    .......
#	    .......
#
#
#    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  W133    W148-3  W158-3  W164    W173    W180-3  W182    W190-3  W191-3  W194    W204-3  W221    W239    W250

package Vcf ;

sub get_line {
    my $self = shift ;
    return $self->{line} ;
}

sub get_INFO {
    my $self = shift ;
    return $self->{INFO} ;
}

sub get_alt_allele_num {
    my $self = shift ;
    my @alleles = split(/\,/,$self->{ALT}) ;

    return scalar @alleles ;
}

sub get_gt_combination_num {
    my $self = shift ;
    my $ploidy = shift ;
    my $alt_allelenum = $self->get_alt_allele_num ;
    my $all_allelenum = $alt_allelenum + 1 ;

    my $output = (($all_allelenum * $all_allelenum - $all_allelenum) / 2 ) + $all_allelenum ;
    $output = 2 if ($ploidy == 1) ;

    return $output ;
}

sub get_b4info {
    my $self = shift ;
    my @list_col2get = qw/CHROM POS ID REF ALT QUAL FILTER/;
    my $lite = {map {$_ => $self->{$_}} @list_col2get };
    bless $lite , "Vcf" ;
    return $lite ;
}

sub get_all {
    my $self = shift ;
    my @list_col2convert = qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT) ;
    push @list_col2convert , @{$self->{samples}} ;

    my $all = {map {$_ => $self->{$_}} @list_col2convert } ;

    my $infotext = $self->{INFO}->{all} ;
    $all->{INFO} = $infotext ;
    $all->{line} = join("\t" , map {$all->{$_}} @list_col2convert ) ;

    bless $all, "Vcf" ;
    return $all ;
}

sub get_basic {
    my $self = shift ;
    my @list_col2convert = qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT) ;

    my $lite = {map {$_ => $self->{$_}} @list_col2convert } ;
    $lite->{line} = join ("\t" , map {$lite->{$_}} @list_col2convert) ;

#    foreach my $col (@list_col2convert) {
#	$lite->{$col} = $self->{$col} ;
#    }

    bless $lite, "Vcf" ;
    return $lite ;
}

1;
