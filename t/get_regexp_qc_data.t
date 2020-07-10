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
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Qc_data}        => [qw{ get_regexp_qc_data }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ get_regexp_qc_data };

diag(   q{Test get_regexp_qc_data from Qc_data.pm v}
      . $MIP::Qc_data::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a valid regexp when file exists
my $data_file_path = catfile( dirname($Bin), q{mip} );
my $recipe_name    = q{mip};
my $regexp = q{perl -nae 'if ($_ =~ /\bour\s\$VERSION\b/xms) {print q{Got version};}'};
my $regexp_key = q{test_add_data_from_regexp};

my @regexp_returns = get_regexp_qc_data(
    {
        data_file_path => $data_file_path,
        regexp         => $regexp,
    }
);

my @expected_qc_data = q{Got version};

## Then return array with single value
is_deeply( \@regexp_returns, \@expected_qc_data,
    q{Got regexp return data from system call} );

## Given a invalid regexp when file exist
my @response = trap { get_regexp_qc_data(
    {
        data_file_path => $data_file_path,
        regexp         => q{perl -e 'print STDERR q{Testing catching STDERR}'},
    }
)
};

is( $response[0], undef, q{Return undef when STDERR} );
like( $trap->stderr, qr/Testing\s+catching\s+STDERR/xms, q{Throw warning when STDERR} );

done_testing();
