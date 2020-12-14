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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Processmanagement::Slurm_processes} => [qw{ submit_slurm_recipe }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Slurm_processes qw{ submit_slurm_recipe };

diag(   q{Test submit_slurm_recipe from Slurm_processes.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

## Given a mock slurm and script
my $case_id = q{case1};
my @dependency_methods =
  qw{ add_to_all case_to_island island island_to_sample island_to_samples sample_to_case sample_to_case_parallel sample_to_island sample_to_sample sample_to_sample_parallel };
my %job_id                       = test_mip_hashes( { mip_hash_name => q{job_id}, } );
my $job_id_chain                 = q{MAIN};
my %max_parallel_processes_count = (
    sample1 => 0,
    sample2 => 0,
);
my @parallel_chains = qw{ OTHER };
my $sample_id       = q{sample1};
my @sample_ids      = qw{ sample1 sample2 };
my $slurm_mock_cmd  = catfile( $Bin, qw{ data modules slurm-mock.pl } );
my $sbatch_file_name =
  catfile( $Bin, qw{ data 643594-miptest test_script fastqc_ADM1059A1.0.sh } );
my $sbatch_script_tracker = 0;

## Given parameters and dependency methods
DEPENDENCY_METHOD:
foreach my $dependency_method (@dependency_methods) {

    my $is_ok = submit_slurm_recipe(
        {
            base_command                      => $slurm_mock_cmd,
            case_id                           => $case_id,
            dependency_method                 => $dependency_method,
            job_id_chain                      => $job_id_chain,
            job_id_href                       => \%job_id,
            log                               => $log,
            max_parallel_processes_count_href => \%max_parallel_processes_count,
            parallel_chains_ref               => \@parallel_chains,
            sample_id                         => $sample_id,
            sample_ids_ref                    => \@sample_ids,
            recipe_file_path                  => $sbatch_file_name,
            recipe_files_tracker              => $sbatch_script_tracker,
        }
    );
## Then recipe shouild be submitted
    ok( $is_ok, q{Submitted recipe } . $dependency_method );
}

done_testing();
