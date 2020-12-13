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
use MIP::Test::Fixtures qw{ test_mip_hashes test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Processmanagement::Slurm_processes} =>
          [qw{ slurm_submit_job_sample_id_dependency_case_dead_end }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Slurm_processes
  qw{ slurm_submit_job_sample_id_dependency_case_dead_end };

diag(
    q{Test slurm_submit_job_sample_id_dependency_case_dead_end from Slurm_processes.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given a mock slurm and script
my $case_id        = q{case1};
my %job_id         = test_mip_hashes( { mip_hash_name => q{job_id}, } );
my $path           = q{MAIN};
my @sample_ids     = qw{ sample1 sample2 };
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );
my $sbatch_file_name =
  catfile( $Bin, qw{ data 643594-miptest test_script fastqc_ADM1059A1.0.sh } );

my %max_parallel_processes_count = (
    sample1 => 1,
    sample2 => 0,
);

## When parallel case chain job id exits
@{ $job_id{case1_MAIN}{case1_parallel_MAIN0} } = qw{ job_id_10 };

slurm_submit_job_sample_id_dependency_case_dead_end(
    {
        base_command                      => $slurm_mock_cmd,
        case_id                           => $case_id,
        job_id_href                       => \%job_id,
        max_parallel_processes_count_href => \%max_parallel_processes_count,
        path                              => $path,
        sample_ids_ref                    => \@sample_ids,
        sbatch_file_name                  => $sbatch_file_name,
    }
);

## Then add job_id returned to ALL
my $expected_return = q{1234};
is( $job_id{ALL}{ALL}[0], $expected_return, q{Added job_id to ALL } );

## Then previous case ids should be cleared and submitted job id added
is( scalar @{ $job_id{case1_MAIN}{case1_MAIN} }, 0, q{Cleared dependencies for case} );

done_testing();
