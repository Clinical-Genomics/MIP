package MIP::Check::Cluster;

#### Copyright 2017 Henrik Stranneheim

use strict;
use warnings;
use warnings qw(FATAL utf8);
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(check_max_core_number);
}

sub check_max_core_number {

##check_max_core_number

##Function : Limit number of cores requested to the maximum number of cores available per node.
##Returns  : "$core_number_requested|$max_cores_per_node"
##Arguments: $max_cores_per_node, $core_number_requested
##         : $max_cores_per_node    => The max number of cores per node
##         : $core_number_requested => The number of cores requested to allocate

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $max_cores_per_node;
    my $core_number_requested;

    my $tmpl = {
        max_cores_per_node => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$max_cores_per_node
        },
        core_number_requested => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$core_number_requested
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    # Check number of cores according to max cores per node
    if ( $core_number_requested > $max_cores_per_node ) {

        # Return maximum available cores
        return $max_cores_per_node;
    }
    else {

        # Core number requested is lower than available cores per node
        return $core_number_requested;
    }
}

1;
