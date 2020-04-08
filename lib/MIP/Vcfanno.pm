package MIP::Vcfanno;

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

## MIPs lib
use MIP::Constants qw{ $LOG_NAME };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_vcfanno_toml parse_toml_config_parameters };
}

sub check_vcfanno_toml {

## Function : Check that the supplied vcfanno toml config has mandatory keys and file exists for annotation array
## Returns  :
## Arguments: $parameter_name    => Name of parameter
##          : $vcfanno_file_toml => Toml config file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_name;
    my $vcfanno_file_toml;

    my $tmpl = {
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        vcfanno_file_toml => {
            defined     => 1,
            required    => 1,
            store       => \$vcfanno_file_toml,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_and_index_existance };
    use MIP::Io::Read qw{ read_from_file };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %vcfanno_config = read_from_file(
        {
            format => q{toml},
            path   => $vcfanno_file_toml,
        }
    );
    my %vcfanno_feature = (
        vcf   => [qw{ file fields ops }],
        other => [qw{ columns file names ops }],
    );

    my $err_msg = q{ is not defined or empty vcfanno toml features. Please check file: }
      . $vcfanno_file_toml;
    my @missing_annotations;

  ANNOTATION:
    foreach my $annotation_href ( @{ $vcfanno_config{annotation} } ) {

        my $file_format = q{other};
        if ( $annotation_href->{file} =~ qr/[.]vcf | [.]vcf[.]gz/xsm ) {
            $file_format = q{vcf};
        }

      FEATURE:
        foreach my $feature ( @{ $vcfanno_feature{$file_format} } ) {

            ## Check mandatory feature keys for vcfanno
            next FEATURE if ( defined $annotation_href->{$feature} );

            push @missing_annotations, q{Feature: } . $feature . $err_msg;
        }
        if (@missing_annotations) {

            ## Broadcast missing files
            foreach my $error_msg (@missing_annotations) {
                $log->fatal($error_msg);
            }
            exit 1;
        }

        ## Check path object exists
        check_filesystem_objects_and_index_existance(
            {
                object_name    => q{file},
                object_type    => q{file},
                parameter_name => $parameter_name,
                path           => $annotation_href->{file},
            }
        );
    }
    return 1;
}

sub parse_toml_config_parameters {

## Function : Parse parameters with TOML config files
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis

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

    ## Check that the supplied vcfanno toml config has mandatory keys and file exists for annotation array
    my %toml_config_parameter = (
        frequency_filter => q{fqa_vcfanno_config},
        sv_annotate      => q{sv_fqa_vcfanno_config},
    );

  CONFIG_FILE:
    while ( my ( $recipe_name, $parameter_name ) = each %toml_config_parameter ) {

        next CONFIG_FILE if ( not $active_parameter_href->{$recipe_name} );

        check_vcfanno_toml(
            {
                parameter_name    => $parameter_name,
                vcfanno_file_toml => $active_parameter_href->{$parameter_name},
            }
        );
    }
    return 1;
}

1;
