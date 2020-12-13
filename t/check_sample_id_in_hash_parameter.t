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
        q{MIP::Active_parameter} => [qw{ check_sample_id_in_hash_parameter }],
        q{MIP::Test::Fixtures}   => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ check_sample_id_in_hash_parameter };

diag(   q{Test check_sample_id_in_hash_parameter from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

my %active_parameter = ( sample_ids => [qw{ sample-1 sample-2 }], );
my %parameter;

## Given no hash parameter
my $is_skipped = check_sample_id_in_hash_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_names_ref   => [qw{ analysis_type }],
        parameter_href        => \%parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);

## Then we have nothing to test
ok( $is_skipped, q{Skipped parameter} );

## Given a hash parameter, when undef sample_ids for parameter
%active_parameter = (
    analysis_type => undef,
    sample_ids    => [qw{ sample-1 sample-2 }],
);

$is_skipped = check_sample_id_in_hash_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_names_ref   => [qw{ analysis_type }],
        parameter_href        => \%parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);

## Then we still have nothing to test
ok( $is_skipped, q{Skipped undef sample id for parameter} );

## Given a hash parameter, when parameter is mandatory
%active_parameter = (
    analysis_type => { q{sample-1} => q{wgs}, },
    sample_ids    => [qw{ sample-1 sample-2 }],
);

trap {
    check_sample_id_in_hash_parameter(
        {
            active_parameter_href => \%active_parameter,
            parameter_names_ref   => [qw{ analysis_type }],
            parameter_href        => \%parameter,
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if mandatory parameter} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message for mandatory parameter} );

## Given a hash parameter, when parameter is not mandatory
%active_parameter = (
    analysis_type => { q{sample-1} => q{wgs}, },
    sample_ids    => [qw{ sample-1 sample-2 }],
);

%parameter = ( analysis_type => { mandatory => q{no}, }, );

$is_skipped = check_sample_id_in_hash_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_names_ref   => [qw{ analysis_type }],
        parameter_href        => \%parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);

## Then we can skip the check
ok( $is_skipped, q{Skipped check for not mandatory parameter} );

## Given a hash parameter, when for ok parameters
%active_parameter = (
    analysis_type => {
        q{sample-1} => q{wgs},
        q{sample-2} => q{wes},
    },
    sample_ids => [qw{ sample-1 sample-2 }],
);

%parameter = ( analysis_type => { mandatory => q{yes}, }, );

my $is_ok = check_sample_id_in_hash_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_names_ref   => [qw{ analysis_type }],
        parameter_href        => \%parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);

## Then we can skip the check
ok( $is_ok, q{Passed check for hash parameter} );
done_testing();
