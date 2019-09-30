package MIP::File::Format::Reference;

use 5.026;
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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ write_references };
}

sub write_references {

## Function : Write references for this analysis to yaml
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $outfile_path          => Outfile path for reference yaml file
##          : $parameter_href        => Holds all parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $outfile_path;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Yaml qw{ write_yaml };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
        }
    );
    my %reference;

  PARAMETER:
    foreach my $parameter_name ( keys %{$active_parameter_href} ) {

        ## Only defined reference parameters
        if ( exists $parameter_href->{$parameter_name}{is_reference} ) {

            if ( ref $active_parameter_href->{$parameter_name} eq q{HASH} ) {
                $reference{$parameter_name} =
                  \%{ $active_parameter_href->{$parameter_name} };
                next PARAMETER;
            }
            if ( ref $active_parameter_href->{$parameter_name} eq q{ARRAY} ) {
                $reference{$parameter_name} =
                  \@{ $active_parameter_href->{$parameter_name} };
                next PARAMETER;
            }

            $reference{$parameter_name} = $active_parameter_href->{$parameter_name};

        }
    }

    # Writes a YAML hash to file
    write_yaml(
        {
            yaml_href      => \%reference,
            yaml_file_path => $outfile_path,
        }
    );
    $log->info( q{Wrote reference YAML file to: } . $outfile_path );

    return;
}

1;
