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

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_toml_annotation_for_tags
      check_vcfanno_toml
      parse_toml_config_parameters
    };
}

sub check_toml_annotation_for_tags {

## Function : Check that TOML annotation contains necessary vcf tags
## Returns  :
## Arguments: $annotation_href      => TOML annotation {REF}
##          : $bcftools_binary_path => Path to bcftools binary
##          : $missing_tag_href     => Files missing vcf tags {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotation_href;
    my $bcftools_binary_path;
    my $missing_tag_href;

    my $tmpl = {
        annotation_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$annotation_href,
            strict_type => 1,
        },
        bcftools_binary_path => {
            defined     => 1,
            required    => 1,
            store       => \$bcftools_binary_path,
            strict_type => 1,
        },
        missing_tag_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$missing_tag_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Vcf qw{ get_vcf_header_line_by_id };

    my $annotation_file_path = $annotation_href->{file};

  VCF_ID_TAG:
    foreach my $vcf_id_tag ( @{ $annotation_href->{fields} } ) {

        my $header_id_line = get_vcf_header_line_by_id(
            {
                bcftools_binary_path => $bcftools_binary_path,
                header_id            => $vcf_id_tag,
                vcf_file_path        => $annotation_file_path,
            }
        );

        next VCF_ID_TAG if defined $header_id_line;

        push @{ $missing_tag_href->{$annotation_file_path} }, $vcf_id_tag;
    }

    return 1;
}

sub check_vcfanno_toml {

## Function : Check that the supplied vcfanno toml config has mandatory keys and file exists for annotation array
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $vcfanno_config_name   => Name of vcfanno config
##          : $vcfanno_functions     => Name of vcfanno functions

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $vcfanno_config_name;
    my $vcfanno_functions;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        vcfanno_config_name => {
            defined     => 1,
            required    => 1,
            store       => \$vcfanno_config_name,
            strict_type => 1,
        },
        vcfanno_functions => {
            defined     => 1,
            required    => 1,
            store       => \$vcfanno_functions,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_and_index_existance };
    use MIP::Io::Read qw{ read_from_file };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $vcfanno_file_toml = $active_parameter_href->{$vcfanno_config_name};
    my %vcfanno_config    = read_from_file(
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
                parameter_name => $vcfanno_config_name,
                path           => $annotation_href->{file},
            }
        );
    }

    ## Check for function file
    if ( $vcfanno_config{functions} and $vcfanno_config{functions}{file} ) {

        ## Check path object exists
        check_filesystem_objects_and_index_existance(
            {
                object_name    => q{file},
                object_type    => q{file},
                parameter_name => $vcfanno_functions,
                path           => $vcfanno_config{functions}{file},
            }
        );
        ## Set path in active parameter
        $active_parameter_href->{$vcfanno_functions} = $vcfanno_config{functions}{file};
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
        variant_annotation => [qw{ vcfanno_config vcfanno_functions }],
        sv_annotate        => [qw{ sv_vcfanno_config sv_vcfanno_functions }],
    );

  RECIPE:
    while ( my ( $recipe_name, $parameter_names_ref ) = each %toml_config_parameter ) {

        next RECIPE if ( not $active_parameter_href->{$recipe_name} );

        my ( $vcfanno_config_name, $vcfanno_functions ) = @{$parameter_names_ref};

        check_vcfanno_toml(
            {
                active_parameter_href => $active_parameter_href,
                vcfanno_config_name   => $vcfanno_config_name,
                vcfanno_functions     => $vcfanno_functions,
            }
        );
    }
    return 1;
}

1;
