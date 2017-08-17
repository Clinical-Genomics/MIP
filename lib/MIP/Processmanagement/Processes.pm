package MIP::Processmanagement::Processes;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    # Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];

use FindBin qw($Bin);    # Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw(unix_standard_streams);
use MIP::Unix::Write_to_file qw(unix_write_to_file);

BEGIN {
    use base qw (Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(print_wait);

}

sub print_wait {

##print_wait

##Function : Calculates when to print "wait" statement and prints "wait" to supplied FILEHANDLE when adequate.
##Returns  : "$process_batches_count"
##Arguments: $process_counter, $max_process_number, $process_batches_count, $FILEHANDLE
##         : $process_counter       => The number of started processes
##         : $max_process_number    => The maximum number of processes to be use before printing "wait" statement
##         : $process_batches_count => Scales the number of $max_process_number processs used after each print "wait" statement
##         : $FILEHANDLE            => FILEHANDLE to print "wait" statment to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $process_counter;
    my $max_process_number;
    my $process_batches_count;
    my $FILEHANDLE;

    my $tmpl = {
        process_counter => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$process_counter
        },
        max_process_number => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$max_process_number
        },
        process_batches_count => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$process_batches_count
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    use MIP::Gnu::Bash qw(gnu_wait);

# Using only nr of processs eq the maximum number of process scaled by the batch count
    if ( $process_counter == $process_batches_count * $max_process_number ) {

        # Print wait statement to filehandle
        gnu_wait( { FILEHANDLE => $FILEHANDLE, } );
        say $FILEHANDLE "\n";

# Increase the maximum number of processs allowed to be used since "wait" was just printed
        $process_batches_count = $process_batches_count + 1;
    }
    return $process_batches_count;
}

1;
