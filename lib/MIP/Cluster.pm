package MIP::Cluster;

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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(get_core_number update_core_number_to_seq_mode);
}

sub get_core_number {

##get_core_number

##Function : Get core number depending on user supplied input exists or not and max number of cores.
##Returns  : "$core_number"
##Arguments: $module_core_number, $modifier_core_number, $max_cores_per_node
##         : $module_core_number   => User input module core numbers to use
##         : $modifier_core_number => Modifier core number dependent on mode of operation of command
##         : $max_cores_per_node   => The max number of cores per node

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $module_core_number;
    my $modifier_core_number;
    my $max_cores_per_node;

    my $tmpl = {
        module_core_number =>
          { strict_type => 1, store => \$module_core_number },
        modifier_core_number => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$modifier_core_number
        },
        max_cores_per_node => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$max_cores_per_node
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    use MIP::Check::Cluster qw(check_max_core_number);

    my $core_number;

    if (   ( defined $module_core_number )
        && ($module_core_number) )
    {

        $core_number = $module_core_number;
    }
    else {

        $core_number = $modifier_core_number;
    }

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            max_cores_per_node    => $max_cores_per_node,
            core_number_requested => $core_number,
        }
    );

    return $core_number;
}

sub update_core_number_to_seq_mode {

##update_core_number_to_seq_mode

##Function : Update the number of cores to be used in the analysis according to sequencing mode requirements.
##Returns  : "$core_number"
##Arguments: $core_number, sequence_run_type
##         : $core_number       => Number of cores to use in the analysis
##         : $sequence_run_type => Type of sequencing [paired-end|single-end]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $core_number;
    my $sequence_run_type;

    my $tmpl = {
        core_number => {
            required    => 1,
            strict_type => 1,
            store       => \$core_number
        },
        sequence_run_type => {
            required    => 1,
            strict_type => 1,
            store       => \$sequence_run_type
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    # Second read direction if present
    if ( $sequence_run_type eq 'paired-end' )
    {

      # 2 processes per file
      $core_number = $core_number + 2;
    }
    elsif ($sequence_run_type eq 'single-end') {

        # Only 1 file and one process
        $core_number = $core_number + 1;
    }
    return $core_number;
}

1;
