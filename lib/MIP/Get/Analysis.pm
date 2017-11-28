package MIP::Get::Analysis;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## Third party module(s)
use autodie;
use List::MoreUtils qw{ all any };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_overall_analysis_type print_program };

}

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};

sub get_overall_analysis_type {

## Function : Detect if all samples has the same sequencing type and return consensus or mixed
## Returns  : q{consensus} | q{mixed} - analysis_type
## Arguments: $analysis_type_href => Analysis_type hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type_href;

    my $tmpl = {
        analysis_type_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$analysis_type_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @analysis_types = (qw{ wes wgs wts });

  ANALYSIS:
    foreach my $analysis_type (@analysis_types) {

        ## If consensus is reached
        if ( all { $_ eq $analysis_type } values %{$analysis_type_href} ) {

            return $analysis_type;
        }
    }

    # No consensus, then it must be mixed
    return q{mixed};
}

sub print_program {

## Function : Print all supported programs in '-ppm' mode
## Returns  :
## Arguments: $define_parameters_file => MIPs define parameters file
##          : $parameter_href         => Parameter hash {REF}
##          : $print_program_mode     => Mode to run modules in

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    ## Default(s)
    my $define_parameters_file;
    my $print_program_mode;

    my $tmpl = {
        define_parameters_file => {
            default =>
              catfile( $Bin, qw{ definitions define_parameters.yaml } ),
            strict_type => 1,
            store       => \$define_parameters_file
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        print_program_mode => {
            default => $arg_href->{print_program_mode} //= 2,
            allow       => [ undef, 0, 1, 2 ],
            strict_type => 1,
            store       => \$print_program_mode
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Yaml qw{ order_parameter_names };
    use MIP::Set::Parameter qw{ set_dynamic_parameter };

    my @printed_programs;

    set_dynamic_parameter(
        {
            parameter_href => $parameter_href,
            aggregates_ref => [q{type:program}],
        }
    );

    ## Adds the order of first level keys from yaml file to array
    my @order_parameters = order_parameter_names(
        {
            file_path => $define_parameters_file,
        }
    );

  PARAMETER:
    foreach my $parameter (@order_parameters) {

        ## Only process programs
        if (
            any { $_ eq $parameter }
            @{ $parameter_href->{dynamic_parameter}{program} }
          )
        {

            if (
                not $parameter =~
                / pbamcalibrationblock | pvariantannotationblock /xsm )
            {

                print {*STDOUT} q{--}
                  . $parameter
                  . $SPACE
                  . $print_program_mode
                  . $SPACE;

                push @printed_programs,
                  $parameter . $SPACE . $print_program_mode;
            }
        }
    }
    print {*STDOUT} $NEWLINE;
    return @printed_programs;
}

1;
