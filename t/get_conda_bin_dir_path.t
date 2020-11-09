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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

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
        q{MIP::Environment::Path} => [qw{ get_conda_bin_dir_path }],
        q{MIP::Test::Fixtures}    => [qw{ test_mip_hashes test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Path qw{ get_conda_bin_dir_path };

diag(   q{Test get_conda_bin_dir_path from Path.pm v}
      . $MIP::Environment::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given a conda bin and bwa bin file path, when exists
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );

my $binary      = q{bwa};
my $recipe_name = q{bwa_mem};
$active_parameter{conda_path} = catfile( $Bin, qw{ data modules miniconda } );

my $is_ok = get_conda_bin_dir_path(
    {
        active_parameter_href => \%active_parameter,
        bin_file              => $binary,
        environment_key       => $recipe_name,
    }
);

ok( $is_ok, q{Found dynamic conda bin path} );

done_testing();
