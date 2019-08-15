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
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Get::Parameter} => [qw{ get_dynamic_conda_path }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_dynamic_conda_path };

diag(   q{Test get_dynamic_conda_path from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a binary and recipe name, but no actual conda bin path
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );

my $binary      = q{bwa};
my $recipe_name = q{bwa_mem};

my $ret = get_dynamic_conda_path(
    {
        active_parameter_href => \%active_parameter,
        bin_file              => $binary,
        conda_bin_file        => q{not_a_conda_bin},
        environment_key       => $recipe_name,
    }
);

## Then exit and throw log message
like(
    $ret,
    qr/Failed \s+ to \s+ find \s+ default \s+ conda \s+ path/xms,
    q{Throw log message if no default conda path}
);

## Given a conda bin and bwa bin file path, when exists
$active_parameter{conda_path} = catfile( $Bin, qw{ data modules miniconda } );

my $bin_path = get_dynamic_conda_path(
    {
        active_parameter_href => \%active_parameter,
        bin_file              => $binary,
        environment_key       => $recipe_name,
    }
);

# Expected path
my $bwa_bin_path = catfile( $active_parameter{conda_path}, qw{ envs test bin } );

is( $bin_path, $bwa_bin_path, q{Found dynamic conda bin path} );

## Given a conda bin when bwa bin file path does not exist
my $faulty_binary = q{not_a_bin_file};

my $err_msg = get_dynamic_conda_path(
    {
        active_parameter_href => \%active_parameter,
        bin_file              => $faulty_binary,
        environment_key       => $recipe_name,
    }
);

like(
    $err_msg,
    qr/Failed \s+ to \s+ find \s+ default \s+ path/xms,
    q{Could not find default binary file}
);

done_testing();
