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
        q{MIP::Set::Parameter} => [qw{ set_recipe_mode }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Parameter qw{ set_recipe_mode };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

diag(   q{Test set_recipe_mode from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given input to set recipe mode to dry_run
Readonly my $TWO => 2;
my $log              = test_log( {} );
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
    }
);
my @recipes = qw{ salmon_quant };

trap {
    set_recipe_mode(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            mode                  => $TWO,
            recipes_ref           => \@recipes,
        }
    )
};

## Then set salmon qunat mode to 2
is( $active_parameter{salmon_quant}, 2, q{Set recipe mode} );
like( $trap->stderr, qr/Set\ssalmon_quant\sto:\s2/xms, q{Write to log} );

done_testing();
