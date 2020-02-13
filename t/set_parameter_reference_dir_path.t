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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ set_parameter_reference_dir_path }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_parameter_reference_dir_path };

diag(   q{Test set_parameter_reference_dir_path from Parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $reference_dir = catdir(qw{ test ref_dir });

my %parameter = (
    file  => q{test},
    array => q{test},
    hash  => q{test},
);

my %active_parameter = (
    reference_dir => $reference_dir,
    file          => q{file_0},
    array         => [qw{ file_1 file_2 }],
    hash          => { file_3 => q{info_key_1} },
);

PARAMETER:
foreach my $parameter_name ( keys %parameter ) {

    set_parameter_reference_dir_path(
        {
            active_parameter_href => \%active_parameter,
            parameter_name        => $parameter_name,
        }
    );
}

## File
is(
    $active_parameter{file},
    catfile( $reference_dir, q{file_0} ),
    q{Set file reference path}
);

is(
    $active_parameter{array}[0],
    catfile( $reference_dir, q{file_1} ),
    q{Set array reference path}
);

UPDATED_FILE:
foreach my $updated_file ( keys %{ $active_parameter{hash} } ) {

    is( $updated_file, catfile( $reference_dir, q{file_3} ), q{Set hash reference path} );

}

done_testing();
