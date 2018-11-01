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
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Get::Parameter} => [qw{ get_recipe_parameters }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_recipe_parameters };

diag(   q{Test get_recipe_parameters from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $recipe_name = q{chanjo_sexcheck};

my %active_parameter = (
    recipe_time        => { chanjo_sexcheck => 1, },
    recipe_core_number => { chanjo_sexcheck => 1, },
    load_env           => {
        mip_env => {
            mip    => undef,
            method => q{conda},
        },
    },
);

my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
    }
);

is( $core_number, 1, q{Got module core_number} );

is( $time, 1, q{Got module time} );

is_deeply(
    \@source_environment_cmds,
    [qw{source activate mip_env}],
    q{Got source environment command}
);

## Test module specific source command
%active_parameter = (
    recipe_time        => { chanjo_sexcheck => 1, },
    recipe_core_number => { chanjo_sexcheck => 1, },
    load_env           => {
        q{py3.6} => {
            chanjo_sexcheck => undef,
            method          => q{conda},
        },
    },
);

( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
    }
);

is_deeply(
    \@source_environment_cmds,
    [qw{ source activate py3.6 }],
    q{Got load env command}
);

done_testing();
