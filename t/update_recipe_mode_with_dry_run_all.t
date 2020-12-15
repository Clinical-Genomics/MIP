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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ update_recipe_mode_with_dry_run_all }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{update_recipe_mode_with_dry_run_all};

diag(   q{Test update_recipe_mode_with_dry_run_all from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

### No update of recipe parameters

my $dry_run_mode = 0;

my @recipes = qw{fastqc_ar bwa_mem peddy_ar};

my %active_parameter = (
    bwa_mem   => 1,
    fastqc_ar => 0,
    peddy_ar  => 2,
);

update_recipe_mode_with_dry_run_all(
    {
        active_parameter_href => \%active_parameter,
        dry_run_all           => $dry_run_mode,
        recipes_ref           => \@recipes,
    }
);

is( $active_parameter{fastqc_ar}, 0, q{No update fastqc_ar} );

is( $active_parameter{bwa_mem}, 1, q{No update bwa_mem} );

is( $active_parameter{peddy_ar}, 2, q{No update peddy_ar} );

### Update of recipe parameters

# Set simulation mode
$dry_run_mode = 1;

update_recipe_mode_with_dry_run_all(
    {
        active_parameter_href => \%active_parameter,
        dry_run_all           => $dry_run_mode,
        recipes_ref           => \@recipes,
    }
);

is( $active_parameter{fastqc_ar}, 0, q{No update fastqc_ar} );

is( $active_parameter{bwa_mem}, 2, q{Update bwa_mem to simulation mode} );

is( $active_parameter{peddy_ar}, 2, q{No update peddy_ar} );

done_testing();
