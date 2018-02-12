package MIP::Cluster;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(get_core_number update_core_number_to_seq_mode);
}

sub get_core_number {

## Function : Get core number depending on user supplied input exists or not and max number of cores.
## Returns  : $core_number
## Arguments: $max_cores_per_node   => The max number of cores per node
##          : $module_core_number   => User input module core numbers to use
##          : $modifier_core_number => Modifier core number dependent on mode of operation of command

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $max_cores_per_node;
    my $module_core_number;
    my $modifier_core_number;

    my $tmpl = {
        max_cores_per_node => {
            defined     => 1,
            required    => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        module_core_number =>
          { store => \$module_core_number, strict_type => 1, },
        modifier_core_number => {
            defined     => 1,
            required    => 1,
            store       => \$modifier_core_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Cluster qw(check_max_core_number);

    my $core_number;

    if ( defined $module_core_number
        && $module_core_number )
    {

        $core_number = $module_core_number;
    }
    else {

        $core_number = $modifier_core_number;
    }

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node    => $max_cores_per_node,
        }
    );

    return $core_number;
}

sub update_core_number_to_seq_mode {

## Function : Update the number of cores to be used in the analysis according to sequencing mode requirements.
## Returns  : $core_number
## Arguments: $core_number       => Number of cores to use in the analysis
##          : $sequence_run_type => Type of sequencing [paired-end|single-end]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $core_number;
    my $sequence_run_type;

    my $tmpl = {
        core_number => {
            required    => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        sequence_run_type => {
            allow       => [qw(paired-end single-end)],
            required    => 1,
            store       => \$sequence_run_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    # Second read direction if present
    if ( $sequence_run_type eq q{paired-end} ) {

        # 2 processes per file
        $core_number = $core_number + 2;
    }
    elsif ( $sequence_run_type eq q{single-end} ) {

        # Only 1 file and one process
        $core_number = $core_number + 1;
    }
    return $core_number;
}

1;
