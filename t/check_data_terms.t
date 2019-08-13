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
use Modern::Perl qw{ 2017 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA %SO_CONSEQUENCE_SEVERITY $SPACE };
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
        q{MIP::Vcfparser}      => [qw{ check_data_terms }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ check_data_terms };

diag(   q{Test check_data_terms from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given a undef consequence term
my $consequence_term = undef;

trap {
    check_data_terms(
        {
            data_href          => \%SO_CONSEQUENCE_SEVERITY,
            term               => $consequence_term,
            data_category_name => q{SO},
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if the term is undef found} );
like( $trap->die, qr/Could\s+not\s+find/xms, q{Throw error if the term is undef found} );

## Given a bad consequence term
$consequence_term = q{not a consequence term};

trap {
    check_data_terms(
        {
            data_href          => \%SO_CONSEQUENCE_SEVERITY,
            term               => $consequence_term,
            data_category_name => q{SO},
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if the term cannot be found} );
like( $trap->die, qr/Could\s+not\s+find/xms, q{Throw error if the term cannot be found} );

## Given a correct consequence term
$consequence_term = q{missense_variant};

my $is_ok = check_data_terms(
    {
        data_href          => \%SO_CONSEQUENCE_SEVERITY,
        term               => $consequence_term,
        data_category_name => q{SO},
    }
);

## Then return true
ok( $is_ok, q{Pass term} );

done_testing();
