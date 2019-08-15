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
        q{MIP::Parse::Parameter} => [qw{ parse_prioritize_variant_callers }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Parameter qw{ parse_prioritize_variant_callers };

diag(   q{Test parse_prioritize_variant_callers from Parameter.pm v}
      . $MIP::Parse::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given no structural active callers, when priority string is ok
my %active_parameter = (
    bcftools_mpileup                       => 1,
    gatk_combinevariants_prioritize_caller => q{gatk,bcftools},
    gatk_variantrecalibration              => 1,
);

my %parameter = (
    cache => {
        structural_variant_callers => [qw{ delly }],
        variant_callers            => [qw{ bcftools_mpileup gatk_variantrecalibration}],
    },
);

my $return = parse_prioritize_variant_callers(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_href        => \%parameter,
    }
);

## Then return undef
is( $return, undef, q{No active structural variant callers} );

## Given structural active callers, when priority string is ok
$active_parameter{sv_svdb_merge_prioritize} = q{delly};
$active_parameter{delly}                    = 1;

my $is_ok = parse_prioritize_variant_callers(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_href        => \%parameter,
    }
);

ok( $is_ok, q{Passed parsing} );

done_testing();
