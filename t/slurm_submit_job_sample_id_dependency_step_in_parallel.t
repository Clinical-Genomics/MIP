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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Processmanagement::Slurm_processes} =>
          [qw{ slurm_submit_job_sample_id_dependency_step_in_parallel }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Slurm_processes
  qw{ slurm_submit_job_sample_id_dependency_step_in_parallel };

diag(
q{Test slurm_submit_job_sample_id_dependency_step_in_parallel from Slurm_processes.pm v}
      . $MIP::Processmanagement::Slurm_processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { no_screen => 1, } );

## Given a mock slurm and script
my $case_id        = q{case1};
my %job_id         = test_mip_hashes( { mip_hash_name => q{job_id}, } );
my $path           = q{MAIN};
my $sample_id      = q{sample1};
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );
my $sbatch_file_name =
  catfile( $Bin, qw{ data 643594-miptest test_script fastqc_ADM1059A1.0.sh } );
my $sbatch_script_tracker = 0;

## When pan key for parallel job exists
@{ $job_id{case1_MAIN}{case1_MAIN_sample1_MAIN} } = qw{ job_id_10 };

slurm_submit_job_sample_id_dependency_step_in_parallel(
    {
        base_command          => $slurm_mock_cmd,
        case_id               => $case_id,
        job_id_href           => \%job_id,
        path                  => $path,
        sample_id             => $sample_id,
        sbatch_file_name      => $sbatch_file_name,
        sbatch_script_tracker => $sbatch_script_tracker,
    }
);

## Then add job_id returned to ALL
my $expected_return = q{1234};
is( $job_id{ALL}{ALL}[0], $expected_return, q{Added job_id to ALL } );

## Then previous sample ids should remain as is
my @sample_job_ids = qw{ job_id_1 job_id_2 };
is_deeply( \@{ $job_id{case1_MAIN}{sample1_MAIN} },
    \@sample_job_ids, q{Added dependencies to sample} );

## Then job ids should be added to the sample parallel step dependency
my @sample_parallel_job_ids = qw{ 1234 };
is_deeply( \@{ $job_id{case1_MAIN}{sample1_parallel_MAIN0} },
    \@sample_parallel_job_ids, q{Added dependencies to sample parallel step} );

done_testing();
