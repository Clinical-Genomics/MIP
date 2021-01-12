package MIP::Set::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      set_programs_for_installation
    };
}

sub set_programs_for_installation {

## Function : Process the lists of programs that has been selected for installation
##          : and update the environment packages
## Returns  :
## Arguments: $active_parameter_href => The entire active parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ array_minus };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check that the options supplied are compatible with each other
    if (    ( scalar @{ $active_parameter_href->{skip_programs} } > 0 )
        and ( scalar @{ $active_parameter_href->{select_programs} } > 0 ) )
    {
        $log->fatal(
            q{"--skip_programs" and "--select_programs" are mutually exclusive command line options}
        );
        exit 1;
    }

    ## Set programs to install depending on pipeline
    if ( scalar @{ $active_parameter_href->{select_programs} } == 0 ) {

      PIPELINE:
        foreach my $pipeline ( @{ $active_parameter_href->{pipelines} } ) {

            push @{ $active_parameter_href->{select_programs} },
              @{ $active_parameter_href->{$pipeline} };
        }
    }

    ## Remove programs that are to be skipped
    delete @{ $active_parameter_href->{container} }{ @{ $active_parameter_href->{skip_programs} } };

    ## Special case for mip_scripts
    @{ $active_parameter_href->{select_programs} } = array_minus(
        @{ $active_parameter_href->{select_programs} },
        @{ $active_parameter_href->{skip_programs} }
    );

    ## Remove all non-selected programs
    if ( scalar @{ $active_parameter_href->{select_programs} } > 0 ) {

        my @non_selects = keys %{ $active_parameter_href->{container} };
        @non_selects =
          array_minus( @non_selects, @{ $active_parameter_href->{select_programs} } );
        delete @{ $active_parameter_href->{container} }{@non_selects};
    }

    return;
}

1;
