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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Analysis}       => [qw{ check_analysis_type_compatibility }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Analysis qw{ check_analysis_type_compatibility };

diag(   q{Test check_analysis_type_compatibility from Analysis.pm v}
      . $MIP::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given compatible analysis type
my $is_ok = check_analysis_type_compatibility(
    {
        analysis_type => q{wgs},
        pipeline      => q{rd_dna},
    }
);

## Then return 1
ok( $is_ok, q{Analysis type and pipeline is compatible} );

## Given a incompatible analysis type
trap {
    check_analysis_type_compatibility(
        {
            analysis_type => q{wts},
            pipeline      => q{rd_dna},
        }
    )
};

## Then croak
ok( $trap->exit, q{Exit if analysis type and pipeline is incompatible} );
like( $trap->stderr, qr/FATAL/xms, q{Throw error} );

done_testing();
