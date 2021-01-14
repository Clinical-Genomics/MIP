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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Build::Rtg_prerequisites} => [qw{ build_rtg_prerequisites }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Build::Rtg_prerequisites qw{ build_rtg_prerequisites };

diag(   q{Test build_rtg_prerequisites from Rtg_prerequisites.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given build parameters
my $parameter_build_name = q{rtg_vcfeval_reference_genome};
my $recipe_name          = q{rtg_vcfeval};
my $slurm_mock_cmd       = catfile( $Bin, qw{ data modules slurm-mock.pl } );

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
$active_parameter{$recipe_name} = 1;

## Unique parameter
$active_parameter{rtg_vcfeval_reference_genome} = q{human_genome.fasta};

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);
my %job_id    = test_mip_hashes( { mip_hash_name => q{job_id}, } );
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);

my %sample_info;

my $is_ok = build_rtg_prerequisites(
    {
        active_parameter_href        => \%active_parameter,
        file_info_href               => \%file_info,
        job_id_href                  => \%job_id,
        log                          => $log,
        parameter_href               => \%parameter,
        recipe_name                  => $recipe_name,
        parameter_build_suffixes_ref => \@{ $file_info{$parameter_build_name} },
        profile_base_command         => $slurm_mock_cmd,
        sample_info_href             => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{ Executed build rtg prerequisites} );

done_testing();
