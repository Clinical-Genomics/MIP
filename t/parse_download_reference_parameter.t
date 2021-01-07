#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Download} => [qw{ parse_download_reference_parameter }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Download qw{ parse_download_reference_parameter };

diag(   q{Test parse_download_reference_parameter from Download.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given versions_ref when array ref and scalar ref
my %reference = (
    ref_1 => [qw{ version_0 version_1 }],
    ref_2 => q{version_2},
);

## When parsing the reference versions
parse_download_reference_parameter( { reference_href => \%reference, } );

my %expected_reference = (
    ref_1 => [qw{ version_0 version_1 }],
    ref_2 => [qw{ version_2 }],
);

## Then all reference versions should be array refs and there should be no scalars left
is_deeply( \%reference, \%expected_reference, q{Parsed versions into array references} );

done_testing();
