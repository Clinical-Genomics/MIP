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
use MIP::Constants qw{ $COLON $DOT $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_definition_file
      get_dependency_tree_from_definition_file
      get_first_level_keys_order_from_definition_file
      get_parameter_definition_file_paths
      get_parameter_from_definition_files
    };
}

sub check_definition_file {

## Function : Parse and check the definition parameters file
## Returns  : %parameter
## Arguments: $define_parameters_path            => File defining the supported parameters
##          : $not_required_definition_file_path => Not required parameter keys
##          : $required_definition_file_path     => Required parameter keys

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $define_parameters_path;
    my $not_required_definition_file_path;
    my $required_definition_file_path;

    my $tmpl = {
        define_parameters_path => {
            defined     => 1,
            required    => 1,
            store       => \$define_parameters_path,
            strict_type => 1,
        },
        not_required_definition_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$not_required_definition_file_path,
            strict_type => 1,
        },
        required_definition_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$required_definition_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parameter qw{ check_parameter_hash };
    use MIP::Io::Read qw{ read_from_file };

    ## Loads MIP definition parameters
    my %parameter = read_from_file(
        {
            format => q{yaml},
            path   => $define_parameters_path,
        }
    );

    ## Load required keys and values for parameters
    my %required = read_from_file(
        {
            format => q{yaml},
            path   => $required_definition_file_path,
        }
    );

    ## Load non required keys and values for parameters
    my %not_required = read_from_file(
        {
            format => q{yaml},
            path   => $not_required_definition_file_path,
        }
    );

    check_parameter_hash(
        {
            file_path         => $define_parameters_path,
            not_required_href => \%not_required,
            parameter_href    => \%parameter,
            required_href     => \%required,
        }
    );
    return %parameter;
}

sub get_dependency_tree_from_definition_file {

## Function : Get dependency tree hash from definition file
## Returns  : %dependency_tree
## Arguments: $level => Level of definition to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $level;

    my $tmpl = {
        level => {
            allow => [
                qw{
                  dragen_rd_dna
                  rd_dna
                  rd_dna_vcf_rerun
                  rd_rna }
            ],
            store       => \$level,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Read qw{ read_from_file };

    my $definition_map_file_path =
      catfile( $Bin, qw{ definitions }, $level . $UNDERSCORE . q{initiation_map.yaml} );
    my %dependency_tree = read_from_file(
        {
            format => q{yaml},
            path   => $definition_map_file_path,
        }
    );

    return %dependency_tree;
}

sub get_parameter_definition_file_paths {

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
                  mip
                  not_required
                  rd_dna
                  rd_dna_vcf_rerun
                  rd_rna
                  required
                  }
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
        mip              => [qw{ mip }],
        not_required     => [qw { not_required }],
        rd_dna           => [qw{ mip analyse rd_dna }],
        rd_dna_vcf_rerun => [qw{ mip analyse rd_dna_vcf_rerun }],
        rd_rna           => [qw{ mip analyse rd_rna }],
        required         => [qw { required }],
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

sub get_first_level_keys_order_from_definition_file {

## Function : Adds the order of first level keys from definition file to array
## Returns  : @order_keys
## Arguments: $file_path => File path to yaml file

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Filehandles
    my $filehandle = IO::Handle->new();

    open $filehandle, q{<}, $file_path
      or croak( q{Cannot open}
          . $DOT
          . $SINGLE_QUOTE
          . $file_path
          . $SINGLE_QUOTE
          . $COLON
          . $SPACE
          . $OS_ERROR
          . $NEWLINE );

    my @order_keys =
      _parse_definition_file_first_level_keys( { filehandle => $filehandle, } );

    close $filehandle;
    return @order_keys;
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
                  mip
                  not_required
                  rd_dna
                  rd_dna_vcf_rerun
                  rd_rna
                  required
                  }
            ],
            store       => \$level,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @definition_file_paths =
      get_parameter_definition_file_paths( { level => $level, } );

    ## Not required parameter definition keys to check
    my $not_required_definition_file_path =
      get_parameter_definition_file_paths( { level => q{not_required}, } );

    ## required parameter definition keys to check
    my $required_definition_file_path =
      get_parameter_definition_file_paths( { level => q{required}, } );

    my %parameter;

  DEFINITION_FILE:
    foreach my $definition_file (@definition_file_paths) {

        %parameter = (
            %parameter,
            check_definition_file(
                {
                    define_parameters_path => $definition_file,
                    not_required_definition_file_path =>
                      $not_required_definition_file_path,
                    required_definition_file_path => $required_definition_file_path,
                }
            ),
        );
    }
    return %parameter;
}

sub _parse_definition_file_first_level_keys {

## Function : Get order of first level keys from definition file
## Returns  : @order_keys
## Arguments: $filehandle => Filehandle to read

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;

    my $tmpl = { filehandle => { store => \$filehandle, }, };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Hold the order of the first level keys from definition file
    my @order_keys;

  LINE:
    while ( my $line = <$filehandle> ) {

        chomp $line;

        ## Next line if header
        next LINE if ( $INPUT_LINE_NUMBER == 1 && $line =~ /\A [-]{3}/sxm );

        ## Next line if commment
        next LINE if ( $line =~ /\A [#]{1}/sxm );

        ## First level key
        my ($key) = $line =~ /\A (\w+):/sxm;

        if ($key) {

            push @order_keys, $key;
            next LINE;
        }
    }
    return @order_keys;
}

1;
