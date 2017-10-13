package MIP::Update::Programs;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ update_program_mode_with_dry_run_all update_program_mode update_prioritize_flag };
}

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub update_program_mode_with_dry_run_all {

## Function : Update program mode depending on dry_run_all flag
## Returns  :
## Arguments: $active_parameter_href, $programs_ref, $dry_run_all
##          : $programs_ref          => Programs in MIP
##          : $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $dry_run_all           => Simulation mode

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
            store       => \$active_parameter_href,
        },
        programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$programs_ref,
        },
        dry_run_all => {
            required    => 1,
            defined     => 1,
            allow       => [ 0, 1, 2 ],
            strict_type => 1,
            store       => \$dry_run_all,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Activated simulation mode
    my $simulation_mode = 2;

    if ($dry_run_all) {

      PROGRAM:
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

sub update_program_mode {

##Function : Update program mode depending on analysis run value as some programs are not applicable for e.g. wes
##Returns  :
##Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $programs_ref            => Programs to update {REF}
##         : $consensus_analysis_type => Consensus analysis_type

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $programs_ref;
    my $consensus_analysis_type;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$programs_ref,
        },
        consensus_analysis_type => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$consensus_analysis_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $consensus_analysis_type ne q{wgs} ) {

        my @warning_msgs;

      PROGRAM:
        foreach my $program ( @{$programs_ref} ) {

            ## Update program mode
            $active_parameter_href->{ q{p} . $program } = 0;

            my $warning_msg =
                q{Turned off: }
              . $program
              . q{ as it is not applicable for }
              . $consensus_analysis_type
              . q{ analysis}
              . $NEWLINE;

            push @warning_msgs, $warning_msg;
        }
        return @warning_msgs;
    }
    return;
}

sub update_prioritize_flag {

##Function : Update prioritize flag depending on analysis run value as some programs are not applicable for e.g. wes
##Returns  :
##Arguments: $programs_ref            => Programs to update {REF}
##         : $prioritize_key          => Prioritize key to update
##         : $consensus_analysis_type => Consensus analysis_type

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $programs_ref;
    my $prioritize_key;
    my $consensus_analysis_type;

    my $tmpl = {
        programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$programs_ref,
        },
        prioritize_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$prioritize_key,
        },
        consensus_analysis_type => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$consensus_analysis_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $consensus_analysis_type ne q{wgs} ) {

        ## Split string into array
        my @callers = split /$COMMA/sxm, $prioritize_key;

      CALLER:
        foreach my $caller ( @{$programs_ref} ) {

            ## Remove all wgs specific callers
            @callers = grep { $caller !~ /^$_/sxm } @callers;
        }
        ## Update sv_svdb_merge_prioritize flag
        $prioritize_key = join q{,}, @callers;
    }
    return $prioritize_key;
}

1;
