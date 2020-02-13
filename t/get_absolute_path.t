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
use Test::Trap;

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
        q{MIP::File::Path}     => [qw{ get_absolute_path }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Path qw{ get_absolute_path };

diag(   q{Test get_absolute_path from File.pm v}
      . $MIP::File::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an existing path
my $existing_path  = catfile( $Bin, qw{ data test_data qc_sample_info.yaml } );
my $parameter_name = q{existing_path};

my $is_ok = get_absolute_path(
    {
        parameter_name => $parameter_name,
        path           => $existing_path,
    }
);

## Then
ok( $is_ok, q{Get absolute path} );

## Given an not existing path
my $not_existing_path = catfile(qw{ data test_data qc_sample_info.yaml });

trap {
    get_absolute_path(
        {
            parameter_name => $parameter_name,
            path           => $not_existing_path,
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if the path cannot be found} );
like( $trap->die, qr/Could \s+ not \s+ find \s+ absolute \s+ path/xms, q{Throw error} );

done_testing();
