package MIP::Qc_data;

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
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ set_qc_data_case_recipe_info set_qc_data_case_recipe_version };
}

sub set_qc_data_case_recipe_info {

## Function : Set case recipe arbitrary info in qc_data hash
## Returns  :
## Arguments: $key => Metafile tag
##          : $qc_data_href => Qc_data hash {REF}
##          : $recipe_name  => Recipe to set attributes in
##          : $value      => Version of program executable

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $key;
    my $qc_data_href;
    my $recipe_name;
    my $value;

    my $tmpl = {
      key => {
        defined     => 1,
        required    => 1,
          strict_type => 1,
          store       => \$key,
      },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        value => {
          defined     => 1,
          required    => 1,
            store       => \$value,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set recipe key value pair for arbitrary info on case level
    $qc_data_href->{recipe}{$recipe_name}{$key} = $value;

return;
}

sub set_qc_data_case_recipe_version {

## Function : Set case recipe version in qc_data hash
## Returns  :
## Arguments: $qc_data_href => Qc_data hash {REF}
##          : $recipe_name  => Recipe to set attributes in
##          : $version      => Version of program executable

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $recipe_name;
    my $version;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        version => {
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $version );

    ## Set recipe version
    $qc_data_href->{recipe}{$recipe_name}{version} = $version;
    return;
}
1;
