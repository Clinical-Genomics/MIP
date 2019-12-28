package MIP::Definition;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };
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
use MIP::Constants qw{ $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_definition_file
      get_definition_file_paths
      get_parameter_from_definition_files
    };
}

sub check_definition_file {

## Function : Parse and check the definition parameters file
## Returns  : %parameter
## Arguments: $define_parameters_path             => File defining the supported parameters
##          : $not_mandatory_definition_file_path => Not mandatory parameter keys
##          : $mandatory_definition_file_path     => Mandatory parameter keys

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $define_parameters_path;
    my $not_mandatory_definition_file_path;
    my $mandatory_definition_file_path;

    my $tmpl = {
        define_parameters_path => {
            defined     => 1,
            required    => 1,
            store       => \$define_parameters_path,
            strict_type => 1,
        },
        not_mandatory_definition_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$not_mandatory_definition_file_path,
            strict_type => 1,
        },
        mandatory_definition_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$mandatory_definition_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parameter qw{  check_parameter_hash };
    use MIP::File::Format::Yaml qw{ load_yaml };

    ## Loads a YAML file into an arbitrary hash and returns it.
    my %parameter = load_yaml( { yaml_file => $define_parameters_path, } );

    ## Load mandatory keys and values for parameters
    my %mandatory = load_yaml(
        {
            yaml_file => $mandatory_definition_file_path,
        }
    );

    ## Load non mandatory keys and values for parameters
    my %not_mandatory = load_yaml(
        {
            yaml_file => $not_mandatory_definition_file_path,
        }
    );

    check_parameter_hash(
        {
            file_path          => $define_parameters_path,
            not_mandatory_href => \%not_mandatory,
            mandatory_href     => \%mandatory,
            parameter_href     => \%parameter,
        }
    );
    return %parameter;
}

sub get_definition_file_paths {

## Function : Get definition file path(s) for definition level
## Returns  : $definition_file_path or @definition_file_paths
## Arguments: $level => Level of definition to inherited from

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $level;

    my $tmpl = {
        level => {
            allow => [
                qw{ analyse
                  dragen_rd_dna
                  download
                  download_rd_dna
                  download_rd_rna
                  install
                  install_rd_dna
                  install_rd_rna
                  mandatory
                  mip
                  not_mandatory
                  rd_dna
                  rd_dna_vcf_rerun
                  rd_rna }
            ],
            store       => \$level,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %definition_map = (
        analyse          => [qw{ mip analyse }],
        dragen_rd_dna    => [qw{ mip analyse dragen_rd_dna }],
        download         => [qw{ mip download }],
        download_rd_dna  => [qw{ mip download download_rd_dna }],
        download_rd_rna  => [qw{ mip download download_rd_rna }],
        install_rd_dna   => [qw{ mip install install_rd_dna }],
        install_rd_rna   => [qw{ mip install install_rd_rna }],
        mandatory        => [qw { mandatory }],
        not_mandatory    => [qw { not_mandatory }],
        mip              => [qw{ mip }],
        rd_dna           => [qw{ mip analyse rd_dna }],
        rd_dna_vcf_rerun => [qw{ mip analyse rd_dna_vcf_rerun }],
        rd_rna           => [qw{ mip analyse rd_rna }],
    );

    my @definition_file_paths;

  LEVEL:
    foreach my $level ( @{ $definition_map{$level} } ) {

        push @definition_file_paths,
          catfile( $Bin, q{definitions}, $level . $UNDERSCORE . q{parameters.yaml} );

    }

    ## Only one definition path - return scalar
    return $definition_file_paths[0] if ( scalar @definition_file_paths == 1 );

    ## Multiple definition paths - return array
    return @definition_file_paths;
}

sub get_parameter_from_definition_files {

## Function : Get parameter hash from definition level
## Returns  : $definition_file_path or @definition_file_paths
## Arguments: $level => Level of definition to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $level;

    my $tmpl = {
        level => {
            allow => [
                qw{ analyse
                  dragen_rd_dna
                  download
                  download_rd_dna
                  download_rd_rna
                  install
                  install_rd_dna
                  install_rd_rna
                  mandatory
                  mip
                  not_mandatory
                  rd_dna
                  rd_dna_vcf_rerun
                  rd_rna }
            ],
            store       => \$level,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @definition_file_paths = get_definition_file_paths( { level => $level, } );

    ## Not mandatory parameter definition keys to check
    my $not_mandatory_definition_file_path =
      get_definition_file_paths( { level => q{not_mandatory}, } );

    ## Mandatory parameter definition keys to check
    my $mandatory_definition_file_path =
      get_definition_file_paths( { level => q{mandatory}, } );

    my %parameter;

  DEFINITION_FILE:
    foreach my $definition_file (@definition_file_paths) {

        %parameter = (
            %parameter,
            check_definition_file(
                {
                    define_parameters_path         => $definition_file,
                    mandatory_definition_file_path => $mandatory_definition_file_path,
                    not_mandatory_definition_file_path =>
                      $not_mandatory_definition_file_path,
                }
            ),
        );
    }

    return %parameter;
}

1;
