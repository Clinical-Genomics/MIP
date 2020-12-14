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
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Validate::Case} => [qw{ check_sample_ids }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Validate::Case qw{ check_sample_ids };

diag(   q{Test check_sample_ids from Case.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given case_id, when no sample_ids
my %active_parameter = ( case_id => q{case-1}, );

trap {
    check_sample_ids(
        {
            case_id        => $active_parameter{case_id},
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if no sample_ids can be found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message no sample_ids can be found} );

## Given case_id, when equals sample_ids
%active_parameter = (
    case_id    => q{sample-1},
    sample_ids => [qw{ sample-1 sample-2 }],
);

trap {
    check_sample_ids(
        {
            case_id        => $active_parameter{case_id},
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    );
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if sample_ids and case_id match} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if sample_ids and case_id match} );

## Given case_id, when duplicate sample_ids
%active_parameter = (
    case_id    => q{case-1},
    sample_ids => [qw{ sample-1 sample-2 sample-1 }],
);

trap {
    check_sample_ids(
        {
            case_id        => $active_parameter{case_id},
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if duplicate samples} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message if duplicate samples} );

## Given case_id, when sample_ids contain "_"
%active_parameter = (
    case_id    => q{case-1},
    sample_ids => [qw{ sample_1 sample-2 sample-3 }],
);

trap {
    check_sample_ids(
        {
            case_id        => $active_parameter{case_id},
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if sample_ids contain "_"} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message if sample_ids contain "_"} );

## Given case_id, when sample_ids are correct
%active_parameter = (
    case_id    => q{case-1},
    sample_ids => [qw{ sample-1 sample-2 sample-3 }],
);

my $is_ok = check_sample_ids(
    {
        case_id        => $active_parameter{case_id},
        sample_ids_ref => \@{ $active_parameter{sample_ids} },
    }
);

## Then  and throw FATAL log message
ok( $is_ok, q{Sample ids are ok} );

done_testing();
