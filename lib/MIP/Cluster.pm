package MIP::Cluster;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX qw{ floor };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $LOG_NAME $SPACE };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_recipe_memory_allocation
      get_core_number
      get_parallel_processes
      update_core_number_to_seq_mode
      update_memory_allocation
    };
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

sub get_parallel_processes {

## Function : Limit number of parallel processes to availble cores and memory
## Returns  : $parallel_processes
## Arguments: $core_number               => Number of cores available
##          : $process_memory_allocation => Memory allocation per process
##          : $recipe_memory_allocation  => Available memory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $core_number;
    my $process_memory_allocation;
    my $recipe_memory_allocation;

    my $tmpl = {
        core_number => {
            defined     => 1,
            required    => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        process_memory_allocation => {
            defined     => 1,
            required    => 1,
            store       => \$process_memory_allocation,
            strict_type => 1,
        },
        recipe_memory_allocation => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_memory_allocation,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $parallel_processes =
      floor( $recipe_memory_allocation / $process_memory_allocation );

    ## Check that the number of processes doesn't exceed the number of available cores
    if ( $parallel_processes > $core_number ) {
        $parallel_processes = $core_number;
    }

    return $parallel_processes;
}

sub update_memory_allocation {

## Function : Calculate recipe memory allocation
## Returns  : $recipe_memory_allocation
## Arguments: $node_ram_memory           => Memory available per node
##          : $parallel_processes        => Number of parallel processes
##          : $process_memory_allocation => Memory requested per process

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $node_ram_memory;
    my $parallel_processes;
    my $process_memory_allocation;

    my $tmpl = {
        node_ram_memory => {
            defined     => 1,
            required    => 1,
            store       => \$node_ram_memory,
            strict_type => 1,
        },
        parallel_processes => {
            defined     => 1,
            required    => 1,
            store       => \$parallel_processes,
            strict_type => 1,
        },
        process_memory_allocation => {
            defined     => 1,
            required    => 1,
            store       => \$process_memory_allocation,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $recipe_memory_allocation = $process_memory_allocation * $parallel_processes;

    ## Check that memory is available
    check_recipe_memory_allocation(
        {
            node_ram_memory          => $node_ram_memory,
            recipe_memory_allocation => $recipe_memory_allocation,
        }
    );

    return $recipe_memory_allocation;

}

sub check_recipe_memory_allocation {

## Function : Check the recipe memory allocation
## Returns  : $recipe_memory_allocation
## Arguments: $node_ram_memory                  => Memory available per node
##          : $memory_allocation_per_process    => Memory requested perl process

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $node_ram_memory;
    my $recipe_memory_allocation;

    my $tmpl = {
        node_ram_memory => {
            defined     => 1,
            required    => 1,
            store       => \$node_ram_memory,
            strict_type => 1,
        },
        recipe_memory_allocation => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_memory_allocation,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Log::MIP_log4perl qw{ retrieve_log };

    ## Retrieves the log instead of importing it
    my $log = retrieve_log( { log_name => $LOG_NAME, } );

    ## Check if memory is available
    if ( $recipe_memory_allocation > $node_ram_memory ) {

        $log->warn( q{Requested memory, }
              . $recipe_memory_allocation
              . q{G, is not available. Allocating }
              . $node_ram_memory
              . q{G} );
        $recipe_memory_allocation = $node_ram_memory;
    }

    return $recipe_memory_allocation;

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
