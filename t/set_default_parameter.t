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
        q{MIP::Active_parameter} => [qw{ set_default_parameter }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_default_parameter };

diag(   q{Test set_default_parameter from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a parameter when scalar
my %active_parameter;
my $scalar_parameter_name    = q{scalar_parameter};
my $scalar_parameter_default = q{a_scalar};
set_default_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => $scalar_parameter_name,
        parameter_default     => $scalar_parameter_default,
    }
);

## Then set default for scalar
is( $active_parameter{$scalar_parameter_name},
    $scalar_parameter_default, q{Set default for scalar parameter} );

## Given a parameter when array
my $array_parameter_name     = q{an_array};
my @array_parameter_defaults = qw{ array parameter };
set_default_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => $array_parameter_name,
        parameter_default     => \@array_parameter_defaults,
    }
);

## Then set default for array
is_deeply( \@{ $active_parameter{$array_parameter_name} },
    \@array_parameter_defaults, q{Set default for array parameter} );

## Given a parameter when hash
my $hash_parameter_name    = q{a_hash};
my %hash_parameter_default = ( key => q{pair}, );
set_default_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => $hash_parameter_name,
        parameter_default     => \%hash_parameter_default,
    }
);

## Then set default for hash
is_deeply( \%{ $active_parameter{$hash_parameter_name} },
    \%hash_parameter_default, q{Set default for hash parameter} );

done_testing();
