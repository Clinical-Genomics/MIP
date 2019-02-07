package MIP::Cluster;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use List::Util qw{ min };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX qw{ floor };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_max_core_number
      get_core_number
      get_memory_constrained_core_number
      update_core_number_to_seq_mode
    };
}

sub check_max_core_number {

## Function : Limit number of cores requested to the maximum number of cores available
## Returns  : $core_number
## Arguments: $core_number_requested => Number of cores requested to allocate
##          : $max_cores_per_node    => Max number of cores per node
##          : $recipe_core_number    => Number of cores allocated by the recipe

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $core_number_requested;
    my $max_cores_per_node;
    my $recipe_core_number;

    my $tmpl = {
        core_number_requested => {
            required    => 1,
            defined     => 1,
            store       => \$core_number_requested,
            strict_type => 1,
        },
        max_cores_per_node => {
            defined     => 1,
            required    => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        recipe_core_number => {
            defined     => 1,
            store       => \$recipe_core_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set default core constraint
    my $core_constraint = $max_cores_per_node;

    ## Update $core_constraint
    if ($recipe_core_number) {

        $core_constraint = min $max_cores_per_node, $recipe_core_number;
    }

    # Check number of cores available
    if ( $core_number_requested > $core_constraint ) {

        # Return maximum available cores
        return $core_constraint;
    }
    else {

        # Core number requested is lower than the number of available cores
        return $core_number_requested;
    }
}

sub get_core_number {

## Function : Get core number depending on user supplied input exists or not and max number of cores.
## Returns  : $core_number
## Arguments: $max_cores_per_node   => The max number of cores per node
##          : $modifier_core_number => Modifier core number dependent on mode of operation of command
##          : $recipe_core_number   => User input module core numbers to use

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $max_cores_per_node;
    my $modifier_core_number;
    my $recipe_core_number;

    my $tmpl = {
        max_cores_per_node => {
            defined     => 1,
            required    => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        modifier_core_number => {
            defined     => 1,
            required    => 1,
            store       => \$modifier_core_number,
            strict_type => 1,
        },
        recipe_core_number => { store => \$recipe_core_number, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $core_number;

    if ( defined $recipe_core_number
        && $recipe_core_number )
    {

        $core_number = $recipe_core_number;
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

sub get_memory_constrained_core_number {

## Function : Limit number of cores requested to the maximum number of cores available given constraints
## Returns  : $core_number
## Arguments: $core_number_requested => Number of cores requested to allocate
##          : $max_cores_per_node    => Max number of cores per node
##          : $memory_allocation     => Memory requested
##          : $node_ram_memory       => Memory available per node
##          : $recipe_core_number    => Number of cores allocated by the recipe

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $max_cores_per_node;
    my $memory_allocation;
    my $node_ram_memory;
    my $recipe_core_number;

    my $tmpl = {
        max_cores_per_node => {
            defined     => 1,
            required    => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        memory_allocation => {
            defined     => 1,
            required    => 1,
            store       => \$memory_allocation,
            strict_type => 1,
        },
        node_ram_memory => {
            defined     => 1,
            required    => 1,
            store       => \$node_ram_memory,
            strict_type => 1,
        },
        recipe_core_number => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_core_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Calculate max usable cores given memory allocation
    my $core_number = floor( ( $recipe_core_number * $node_ram_memory ) /
          ( $max_cores_per_node * $memory_allocation ) );

    if ( $core_number < 1 ) {
        croak q{Requested memory not available! Try to increase recipe_core_number};
    }

    ## Limit number of cores to the maximum number of cores available
    $core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node    => $max_cores_per_node,
            recipe_core_number    => $recipe_core_number,
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
            allow       => [qw{ paired-end single-end }],
            required    => 1,
            store       => \$sequence_run_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

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
