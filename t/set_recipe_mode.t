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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ set_recipe_mode }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_recipe_mode };
use MIP::Test::Fixtures qw{ test_log };

diag(   q{Test set_recipe_mode from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given input to set recipe mode to dry_run
my %active_parameter;
my @recipes = qw{ salmon_quant };

trap {
    set_recipe_mode(
        {
            active_parameter_href => \%active_parameter,
            mode                  => 2,
            recipes_ref           => \@recipes,
        }
    )
};

## Then set salmon quant recipe to mode equals 2
is( $active_parameter{salmon_quant}, 2, q{Set recipe mode} );
like( $trap->stderr, qr/Set \s+ salmon_quant \s+ to: \s+ 2/xms, q{Write to log} );

done_testing();
