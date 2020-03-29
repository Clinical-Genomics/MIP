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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Active_parameter} => [qw{ get_active_parameter_attribute }],
        q{MIP::Test::Fixtures}   => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ get_active_parameter_attribute };

diag(   q{Test get_active_parameter_attribute from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a parameter when attribute is uninitilized
my $parameter_name   = q{analysis_type};
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
    }
);
my $sample_id = $active_parameter{sample_ids}[0];

## When no attribute is supplied
my %attribute = get_active_parameter_attribute(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => $parameter_name,
    }
);

## Then return entire hash
is_deeply( \%{ $active_parameter{$parameter_name} },
    \%attribute, q{Got entire attribute hash} );

## When attribute is uninitilized
my $is_ok = get_active_parameter_attribute(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => $parameter_name,
        attribute             => q{is_not initilized},
    }
);

## Then return undef
is( $is_ok, undef, q{Skipped attribute} );

## When hash of scalar
my $is_scalar = get_active_parameter_attribute(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => $parameter_name,
        attribute             => $sample_id,
    }
);

## Then return scalar
ok( $is_scalar, q{Got scalar attribute} );

## When hash of array
$active_parameter{slurm}{email_types} = [qw{ begin fail }];
my @email_types = get_active_parameter_attribute(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => q{slurm},
        attribute             => q{email_types}
    }
);

## Then return array
is_deeply( \@{ $active_parameter{slurm}{email_types} },
    \@email_types, q{Got array attribute} );

## When hash of hash
$active_parameter{slurm}{cmd}{sacct} = { description => q{a slurm command}, };
my %sacct_attribute = get_active_parameter_attribute(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => q{slurm},
        attribute             => q{cmd},
    }
);

## Then return hash
is_deeply( \%{ $active_parameter{slurm}{cmd} }, \%sacct_attribute,
    q{Got hash attribute} );

done_testing();
