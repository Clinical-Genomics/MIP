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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Check::Parameter} => [qw{ check_prioritize_variant_callers }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Parameter qw{ check_prioritize_variant_callers };

diag(   q{Test check_prioritize_variant_callers from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given active callers, when priority string is ok
my %active_parameter = (
    bcftools_mpileup                       => 1,
    gatk_combinevariants_prioritize_caller => q{gatk,bcftools},
    gatk_variantrecalibration              => 1,
);

my %parameter = (
    cache => {
        variant_callers => [qw{ bcftools_mpileup gatk_variantrecalibration}],
    },
);

my $is_ok = check_prioritize_variant_callers(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_href        => \%parameter,
        parameter_name        => q{gatk_combinevariants_prioritize_caller},
        variant_callers_ref   => \@{ $parameter{cache}{variant_callers} },
    }
);

## Then the sub should return true
ok( $is_ok, q{Passed} );

## Given a priority string with missing variant caller
$active_parameter{gatk_combinevariants_prioritize_caller} = q{gatk};

trap {
    check_prioritize_variant_callers(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
            parameter_name        => q{gatk_combinevariants_prioritize_caller},
            variant_callers_ref   => \@{ $parameter{cache}{variant_callers} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if parameter does not contain active variant caller} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if parameter does not contain active variant caller} );

## Given an not activated variant caller
$active_parameter{bcftools_mpileup}                       = 0;
$active_parameter{gatk_combinevariants_prioritize_caller} = q{gatk,bcftools};

trap {
    check_prioritize_variant_callers(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
            parameter_name        => q{gatk_combinevariants_prioritize_caller},
            variant_callers_ref   => \@{ $parameter{cache}{variant_callers} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if parameter does not contains deactive variant caller} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if parameter does not contains deactive variant caller} );

## Given an other variant caller, when not part of priority string
$active_parameter{bcftools_mpileup} = 1;
$active_parameter{gatk_combinevariants_prioritize_caller} =
  q{gatk,bcftools, NOT_A_VARIANT_CALLER};

trap {
    check_prioritize_variant_callers(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
            parameter_name        => q{gatk_combinevariants_prioritize_caller},
            variant_callers_ref   => \@{ $parameter{cache}{variant_callers} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if parameter does not match any supported variant caller} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if does not match any supported variant caller} );

done_testing();
