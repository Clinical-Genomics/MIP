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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Recipes::Build::Bwa_prerequisites} => [qw{ build_bwa_mem2_prerequisites }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Build::Bwa_prerequisites qw{ build_bwa_mem2_prerequisites };

diag(   q{Test build_bwa_mem2_prerequisites from Bwa_prerequisites.pm v}
      . $MIP::Recipes::Build::Bwa_prerequisites::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given build parameters
my $parameter_build_name = q{bwa_build_reference};
my $recipe_name          = q{bwa_mem};
my $slurm_mock_cmd       = catfile( $Bin, qw{ data modules slurm-mock.pl } );

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
## Submission via slurm_mock
$active_parameter{$recipe_name} = 1;

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);
my %job_id;
my %parameter = test_mip_hashes( { mip_hash_name => q{recipe_parameter}, } );

my %sample_info;

my $is_ok = build_bwa_mem2_prerequisites(
    {
        active_parameter_href        => \%active_parameter,
        file_info_href               => \%file_info,
        job_id_href                  => \%job_id,
        log                          => $log,
        parameter_build_suffixes_ref => \@{ $file_info{$parameter_build_name} },
        parameter_href               => \%parameter,
        profile_base_command         => $slurm_mock_cmd,
        recipe_name                  => $recipe_name,
        sample_info_href             => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{ Executed build bwa_mem prerequisites} );

done_testing();
