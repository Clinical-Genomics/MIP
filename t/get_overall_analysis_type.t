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
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Get::Analysis}  => [qw{ get_overall_analysis_type }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Analysis qw{ get_overall_analysis_type };

diag(   q{Test get_overall_analysis_type from Get.pm v}
      . $MIP::Get::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Base arguments
my %analysis_type_consensus = (
    sample_1 => q{wgs},
    sample_2 => q{wgs},
    sample_3 => q{wgs},
);

my %analysis_type_mixed = (
    sample_1 => q{wgs},
    sample_2 => q{wes},
    sample_3 => q{wes},
);

my %analysis_type_faulty = (
    sample_1 => q{wgs},
    sample_2 => q{wes},
    sample_3 => q{not_a_analysis_type},
);

## Given a consensus analysis type
my $consensus_analysis_type = get_overall_analysis_type(
    {
        analysis_type_href => \%analysis_type_consensus,
        log                => $log,
    }
);

## Then a consensus analysis type should be set
is( $consensus_analysis_type, q{wgs}, q{Correct return of consensus analysis typ} );

## Given a mixed analysis type
my $mixed_analysis_type = get_overall_analysis_type(
    {
        analysis_type_href => \%analysis_type_mixed,
        log                => $log,
    }
);

## Then a mixed analysis type should be set
is( $mixed_analysis_type, q{mixed}, q{Correct return of mixed analysis type} );

## Given a faulty analysis type
trap {
    get_overall_analysis_type(
        {
            analysis_type_href => \%analysis_type_faulty,
            log                => $log,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if analysis type is not supported} );
like(
    $trap->stderr,
    qr/is \s+ not \s+ a \s+ supported \s+ analysis_type/xms,
    q{Throw fatal log message if analysis type is not supported}
);

done_testing();
