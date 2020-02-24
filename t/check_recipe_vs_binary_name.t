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
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Parameter} => [qw{ check_recipe_vs_binary_name }],
        q{MIP::Test::Fixtures}   => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ check_recipe_vs_binary_name };

diag(   q{Test check_recipe_vs_binary_name from Parameter.pm v}
      . $MIP::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a recipe with no identical name as a program binary
my %parameter = test_mip_hashes( { mip_hash_name => q{define_parameter}, } );

my @recipe_names = qw{ bwa_mem sv_reformat };

## Check that recipe name and program name are not identical
my $is_ok = check_recipe_vs_binary_name(
    {
        parameter_href   => \%parameter,
        recipe_names_ref => \@recipe_names,
    }
);

## Then return true
ok( $is_ok, q{No identical recipes and program names} );

## Given a recipe with identical name as a program binary
push @{ $parameter{fastqc}{program_executables} }, q{fastqc};

push @recipe_names, q{fastqc};

## Check that recipe name and program name are not identical
trap {
    check_recipe_vs_binary_name(
        {
            parameter_href   => \%parameter,
            recipe_names_ref => \@recipe_names,
        }
    );
};

## Then croak and die
is( $trap->leaveby, q{die}, q{Exit if the recipe and program name are identical} );
like( $trap->die, qr{\AIdentical\s+names\s+for\s+recipe\s+and\s+program}xms,
    q{Throw error} );

done_testing();
