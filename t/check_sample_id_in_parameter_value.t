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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Check::Parameter} => [qw{ check_sample_id_in_parameter_value }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Parameter qw{ check_sample_id_in_parameter_value };

diag(   q{Test check_sample_id_in_parameter_value from Pedigree.pm v}
      . $MIP::Check::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create log
my $log = test_log( {} );

## Given all samples in pedigree present in analysis
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
my %parameter        = test_mip_hashes( { mip_hash_name => q{recipe_parameter}, } );

my $is_skipped = check_sample_id_in_parameter_value(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_names_ref   => [qw{ subject_id }],
        parameter_href        => \%parameter,
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);

## Then return true
ok( $is_skipped, q{Skipped parameter} );

## Given a non-matching sample
$active_parameter{subject_id}{ADM1059A1} = q{ADM1059A1};
$active_parameter{subject_id}{ADM1059A2} = q{not_a_sample};
$active_parameter{subject_id}{ADM1059A3} = q{ADM1059A3};

trap {
    check_sample_id_in_parameter_value(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_names_ref   => [qw{ subject_id }],
            parameter_href        => \%parameter,
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        }
    );
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the sample id has no match } );
like(
    $trap->stderr,
    qr/Could \s+ not \s+ find \s+ matching \s+ sample_id/xms,
    q{Throw fatal log message if the sample id has no match }
);

## Given a missing sample
$active_parameter{subject_id}{ADM1059A2} = undef;

trap {
    check_sample_id_in_parameter_value(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_names_ref   => [qw{ subject_id }],
            parameter_href        => \%parameter,
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        }
    );
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the sample id cannot be found} );
like(
    $trap->stderr,
    qr/Could \s+ not \s+ find \s+ value/xms,
    q{Throw fatal log message if the sample id cannot be found }
);

done_testing();
