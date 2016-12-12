package VarSelect::Log ;

use strict ;
use Carp ;

sub new {
    my $invocant = shift ;
    my $class    = ref $invocant || $invocant ;
    my %initial_attritubes = @_ ;

    my $self = {
        file    => undef ,
        %initial_attritubes ,
    } ;

    open (my $TGT,">$self->{file}") || confess "Error! Can NOT open $self->{file}!" ;

    $self->{filehandle} = $TGT ;

    bless ($self,$class) ;
    return $self ;
}

sub close {
    my $self = shift ;
    close $self->{filehandle} ;
}

sub write { 
    my $self = shift ;
    my $msg = shift ;
    logit(timelog("$msg") , $self->{filehandle} ) ;
    
}

sub andRun { 
    my $self = shift ;
    my $cmd = shift ;
    $self->write($cmd) ;
    `$cmd` ;
}

sub loadlog {
    my $self = shift ;
    my $file_log = shift ;

    open (my $SRC, $file_log) ;
    while (<$SRC>) {
	chomp ;
	logit($_ , $self->{filehandle} ) ;
    }
    close $SRC ;
}

sub timelog {
    my $text = shift ;
    my $time = localtime(time) ;

    return "$time\t$text" ;
}

sub logit {
    my $msg = shift ;
    my $fh = shift ;

    print $fh "$msg\n" ;
}


#sub tlog {  # timelog + logit , use write
#    my $msg = shift ;
#    logit(timelog("$msg"),$LOG) ;
#}

#sub logNrun { # andRun 
#    my $cmd = shift ;
#
#    tlog($cmd) ;
#    `$cmd` ;
#}

1;
