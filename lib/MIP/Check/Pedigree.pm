package MIP::Check::Pedigree;

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
use List::MoreUtils qw { any };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_pedigree_sample_allowed_values };
}

## Constants
Readonly my $SPACE => q{ };

sub check_pedigree_sample_allowed_values {

## Function : Check that the pedigree sample key values are allowed
## Returns  :
## Arguments: $log           => Log object
##          : $pedigree_href => YAML pedigree info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $pedigree_href;

    my $tmpl = {
        log           => { store => \$log },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    #    my $log = Log::Log4perl->get_logger(q{MIP});

    my %allowed_values = (
        analysis_type => [qw{ wes wgs wts cancer}],
        phenotype     => [qw{ affected unaffected unknown }],
        sample_origin => [qw{ normal tumor }],
        sex           => [qw{ female male other unknown }],
    );

    ## Check input to sample_info hash for at sample level
  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

      SAMPLE_KEY:
        foreach my $key ( keys %{$pedigree_sample_href} ) {

            ## No defined allowed values
            next SAMPLE_KEY if ( not exists $allowed_values{$key} );

            ## If element is not part of array
            next SAMPLE_KEY
              if (
                any { $_ eq $pedigree_sample_href->{$key} }
                @{ $allowed_values{$key} }
              );

            $log->fatal(
                q{Pedigree file key: }
                  . $key
                  . q{ found illegal value: }
                  . $pedigree_sample_href->{$key}
                  . q{ allowed values are '}
                  . join q{' '},
                @{ $allowed_values{$key} }
            );
            $log->fatal(q{Please correct the entry before analysis.});
            $log->fatal(q{MIP: Aborting run.});
            exit 1;
        }
    }
    return 1;
}

1;
