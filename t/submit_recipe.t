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
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Processmanagement::Processes} => [qw{ submit_recipe }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ submit_recipe };

diag(   q{Test submit_recipe from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Base arguments
my $case_id = q{case1};
my $path    = q{MAIN};
my $recipe_file_path =
  catfile( $Bin, qw{ data 643594-miptest test_script fastqc_ADM1059A1.0.sh } );
my $sample_id      = q{sample1};
my @sample_ids     = ($sample_id);
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );
my %job_id         = test_mip_hashes( { mip_hash_name => q{job_id}, } );

## Given a sbatch recipe and parameter, when using slurm-mock

my $is_ok = submit_recipe(
    {
        base_command       => $slurm_mock_cmd,
        case_id            => $case_id,
        dependency_method  => q{island_to_samples},
        job_id_href        => \%job_id,
        log                => $log,
        job_id_chain       => $path,
        recipe_file_path   => $recipe_file_path,
        sample_ids_ref     => \@sample_ids,
        submission_profile => q{slurm},
    }
);

## Then return true
ok( $is_ok, q{Submited recipe} );

## Given a sbatch recipe and parameter, when no SLURM installed and faulty project id in sbatch script
trap {
    submit_recipe(
        {
            case_id            => $case_id,
            dependency_method  => q{island_to_samples},
            job_id_href        => \%job_id,
            log                => $log,
            job_id_chain       => $path,
            recipe_file_path   => $recipe_file_path,
            sample_ids_ref     => \@sample_ids,
            submission_profile => q{slurm},
        }
    )
};

## Then
is( $trap->leaveby, q(die), q{Exit if no SLURM or faulty project id} );

done_testing();
