#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Parameter}      => [qw{ parse_reference_path }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ parse_reference_path };

diag(   q{Test parse_reference_path from Parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $reference_dir = catdir(qw{ test ref_dir });

## Given reference parameters when key "reference" exists or not
my %parameter = (
    array => {
        associated_recipe => [qw{ recipe_1 }],
        reference         => 1,
    },
    file_name_1 => {
        associated_recipe => [qw{ recipe_1 }],
        reference         => 1,
    },
    file_name_2 => { associated_recipe => [qw{ recipe_2 }], },
    hash        => {
        associated_recipe => [qw{ recipe_1 }],
        reference         => 1,
    },
);

my %active_parameter = (
    array         => [qw{ file_1 file_2 }],
    file_name_1   => q{file_0},
    file_name_2   => q{file_2},
    hash          => { file_3 => q{info_key_1} },
    recipe_1      => 1,
    reference_dir => $reference_dir,
);

parse_reference_path(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    }
);

## Then do NOT set mip reference dir path for parameter when "reference" key does NOT exists
is( $active_parameter{file_name_2}, q{file_2}, q{Skipped setting file reference path} );

## Then set set mip reference dir path for parameter when "reference" key exists
is(
    $active_parameter{file_name_1},
    catfile( $reference_dir, q{file_0} ),
    q{Set file reference path}
);

## Then set set mip reference dir path for parameter when "reference" key exists
is(
    $active_parameter{array}[0],
    catfile( $reference_dir, q{file_1} ),
    q{Set array reference path}
);
## Then set set mip reference dir path for parameter when "reference" key exists
UPDATED_FILE:

foreach my $updated_file ( keys %{ $active_parameter{hash} } ) {

    is( $updated_file, catfile( $reference_dir, q{file_3} ), q{Set hash reference path} );

}

done_testing();
