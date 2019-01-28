#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = '1.0.0';

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
        q{MIP::Parse::Parameter} => [qw{ parse_start_with_recipe }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Parameter qw{ parse_start_with_recipe };

diag(   q{Test parse_start_with_recipe from Parameter.pm v}
      . $MIP::Parse::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given no defined start_with_recipe parameter
my %active_parameter;
my %parameter = ( bwa_mem => { default => 0, }, );
my $initiation_file =
  catfile( dirname($Bin), qw{ definitions rd_dna_initiation_map.yaml } );

my $return = parse_start_with_recipe(
    {
        active_parameter_href => \%active_parameter,
        initiation_file       => $initiation_file,
        log                   => $log,
        parameter_href        => \%parameter,
    },
);

## Then skip parsing
is( $return, undef, q{Skip parsing} );

## Given start_with_recipe parameter, when defined
$active_parameter{start_with_recipe} = q{bwa_mem};

my $is_ok = parse_start_with_recipe(
    {
        active_parameter_href => \%active_parameter,
        initiation_file       => $initiation_file,
        log                   => $log,
        parameter_href        => \%parameter,
    },
);

## Then return true for successful parsing
ok( $is_ok, q{Parsed programs from start_with_flag} );

done_testing();

