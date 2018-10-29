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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Check::Parameter} => [qw{ check_recipe_exists_in_hash }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Parameter qw{ check_recipe_exists_in_hash };

diag(   q{Test check_recipe_exists_in_hash from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log();

## Given recipe names
my %parameter = ( q{bcftools_mpileup} => 1, );

my %active_parameter = (
    module_time => {
        bwa_mem          => 1,
        bcftools_mpileup => 1,
    },
    associated_recipe => [ qw{ fastqc_ar }, ],
);
## When one does not exist in truth hash
trap {
    check_recipe_exists_in_hash(
        {
            log            => $log,
            parameter_name => q{module_time},
            query_ref      => \%{ $active_parameter{module_time} },
            truth_href     => \%parameter,
        }
      )
};

## Then exist and throw error
ok( $trap->exit, q{Exit if recipe key does not exist} );
like( $trap->stderr, qr/FATAL/xms, q{Throw FATAL log message} );

## Given recipe names
%parameter = (
    q{bcftools_mpileup} => 1,
    q{bwa_mem}          => 1,
);

## When all exists in truth hash
my $return = check_recipe_exists_in_hash(
    {
        log            => $log,
        parameter_name => q{module_time},
        query_ref      => \%{ $active_parameter{module_time} },
        truth_href     => \%parameter,
    }
);
is( $return, undef, q{All recipe keys exists in truth hash} );

## Given recipe names, when none exists in truth hash
trap {
    check_recipe_exists_in_hash(
        {
            log            => $log,
            parameter_name => q{associated_recipe},
            query_ref      => \@{ $active_parameter{associated_recipe} },
            truth_href     => \%parameter,
        }
      )
};

## Then exist and throw error
ok( $trap->exit, q{Exit if recipe element does not exist} );
like( $trap->stderr, qr/FATAL/xms, q{Throw FATAL log message} );

## Given recipe names
%parameter = (
    q{bcftools_mpileup} => 1,
    q{bwa_mem}          => 1,
    q{fastqc_ar}        => 1,
);

## When all exists in truth hash
$return = check_recipe_exists_in_hash(
    {
        log            => $log,
        parameter_name => q{associated_recipe},
        query_ref      => \@{ $active_parameter{associated_recipe} },
        truth_href     => \%parameter,
    }
);
is( $return, undef, q{All recipe element exists in truth hash} );

## Given a recipe name
my $recipe_name = q{bwa_memA};

%parameter = ( bwa_mem => 1, );

## When recipe does not exists in truth hash
trap {
    check_recipe_exists_in_hash(
        {
            log            => $log,
            parameter_name => $recipe_name,
            query_ref      => \$recipe_name,
            truth_href     => \%parameter,
        }
      )
};

## Then exist and throw error
ok( $trap->exit, q{Exit if recipe does not exist} );
like( $trap->stderr, qr/FATAL/xms, q{Throw FATAL log message} );

## Given a recipe name
$recipe_name = q{bwa_mem};

## When recipe exists in truth hash
$return = check_recipe_exists_in_hash(
    {
        log            => $log,
        parameter_name => $recipe_name,
        query_ref      => \$recipe_name,
        truth_href     => \%parameter,
    }
);
is( $return, undef, q{Recipe exists in truth hash} );

done_testing();
