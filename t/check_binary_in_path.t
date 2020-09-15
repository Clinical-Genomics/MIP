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
use MIP::Test::Fixtures qw{ test_log test_standard_cli test_mip_hashes };

my $VERBOSE = 1;
our $VERSION = 1.08;

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
        q{MIP::Environment::Path} => [qw{ check_binary_in_path }],
        q{MIP::Test::Fixtures}    => [qw{ test_log test_standard_cli test_mip_hashes}],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Path qw{ check_binary_in_path };

diag(   q{Test check_binary_in_path from Path.pm v}
      . $MIP::Environment::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given existing binary
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
$active_parameter{conda_path} =
  catdir( dirname($Bin), qw{ t data modules miniconda } );
my $binary       = q{samtools};
my $program_name = q{samtools};

trap {
    check_binary_in_path(
        {
            active_parameter_href => \%active_parameter,
            binary                => $binary,
            program_name          => $program_name,
        }
    )
};

## Then return true
ok( $trap->return, q{Binary is found} );

## Given no existing binary
my $no_binary = q{Nothing to see here};

trap {
    check_binary_in_path(
        {
            active_parameter_href => \%active_parameter,
            binary                => $no_binary,
            program_name          => $program_name,
        }
    )
};

## Then exit and die
is( $trap->leaveby, q{die}, q{Exit if binary cannot be found} );

done_testing();
