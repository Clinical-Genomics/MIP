package MIP::Update::Programs;

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
    our @EXPORT_OK = qw(update_program_mode_with_dry_run_all);
}

sub update_program_mode_with_dry_run_all {

##update_program_mode_with_dry_run_all

##Function : Update program mode depending on dry_run_all flag
##Returns  : ""
##Arguments: $active_parameter_href, $programs_ref, $dry_run_all
##         : $programs_ref          => Programs in MIP
##         : $active_parameter_href => The active parameters for this analysis hash {REF}
##         : $dry_run_all           => Simulation mode

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $programs_ref;
    my $active_parameter_href;
    my $dry_run_all;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$programs_ref
        },
        dry_run_all => {
            required    => 1,
            defined     => 1,
            allow       => [ 0, 1, 2 ],
            strict_type => 1,
            store       => \$dry_run_all
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Activated simulation mode
    my $simulation_mode = 2;

    if ($dry_run_all) {

      PROGRAMS:
        foreach my $program_name ( @{$programs_ref} ) {

            ## If program is activated
            if ( $active_parameter_href->{$program_name} ) {

                # Change program mode to simulation
                $active_parameter_href->{$program_name} = $simulation_mode;
            }
        }
    }
    return;
}

1;
