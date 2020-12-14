#! /usr/bin/env perl

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
        q{MIP::Get::File} => [qw{ get_path_entries }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::File qw{ get_path_entries };

diag(   q{Test get_path_entries from File.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $CORRECT_PATHS => 3;

my @paths;
my %sample_info = (
    first_level_key => {
        path         => catfile(qw{a test file}),
        outdirectory => catfile(qw{a test}),
        outfile      => q{file},
    },
    first_level_key_2 => {
        outdirectory => catfile(qw{a test}),
        outfile      => q{file_2},
    },
    program => { second_level_key   => { path => catfile(qw{a test file_3}) }, },
    program => { second_level_key_2 => { path => catfile(qw{a test file_3}) }, },
);

## Collects all programs file path(s) created by MIP located in %sample_info
get_path_entries(
    {
        sample_info_href => \%sample_info,
        paths_ref        => \@paths,
    }
);
is( scalar @paths, $CORRECT_PATHS, q{Found correct number of paths} );

isnt( $paths[0], undef, q{Found first path} );

isnt( $paths[1], undef, q{Found second path} );

is( $paths[$CORRECT_PATHS], undef, q{No duplicated paths} );

done_testing();
