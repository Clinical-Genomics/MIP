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
        q{MIP::Parameter}      => [qw{ get_cache }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ get_cache };

diag(   q{Test get_cache from Parameter.pm v}
      . $MIP::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a parameter hash
my %parameter = (
    cache => {
        scalar => 1,
        array  => [qw{ red blue}],
        hash   => {
            program  => q{MIP},
            language => q{perl},
        }
    },
);

## When parameter is uninitilized
my $is_ok = get_cache(
    {
        parameter_href => \%parameter,
        parameter_name => q{not_defined},
    }
);

## Then return undef
is( $is_ok, undef, q{Skipped parameter} );

## When scalar
my $is_scalar = get_cache(
    {
        parameter_href => \%parameter,
        parameter_name => q{scalar},
    }
);

## Then return scalar
ok( $is_scalar, q{Got scalar cache} );

## When array
my @colors = get_cache(
    {
        parameter_href => \%parameter,
        parameter_name => q{array},
    }
);

## Then return array
is_deeply( \@{ $parameter{cache}{array} }, \@colors, q{Got array cache} );

## When hash
my %program = get_cache(
    {
        parameter_href => \%parameter,
        parameter_name => q{hash},
    }
);

## Then return hash
is_deeply( \%{ $parameter{cache}{hash} }, \%program, q{Got hash cache} );

done_testing();
