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
          [qw{ slurm_submit_job_no_dependency_add_to_sample }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Slurm_processes
  qw{ slurm_submit_job_no_dependency_add_to_sample };

diag(   q{Test slurm_submit_job_no_dependency_add_to_sample from Slurm_processes.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a mock slurm and script
my $case_id        = q{case1};
my %job_id         = test_mip_hashes( { mip_hash_name => q{job_id}, } );
my $path           = q{MAIN};
my $sample_id      = q{sample1};
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );
my $sbatch_file_name =
  catfile( $Bin, qw{ data 643594-miptest test_script fastqc_ADM1059A1.0.sh } );
my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

slurm_submit_job_no_dependency_add_to_sample(
    {
        base_command     => $slurm_mock_cmd,
        case_id          => $case_id,
        job_id_href      => \%job_id,
        log              => $log,
        path             => $path,
        sample_id        => $sample_id,
        sbatch_file_name => $sbatch_file_name,
    }
);

## Then add job_id returned to ALL
my $expected_return = q{1234};
is( $job_id{ALL}{ALL}[0], $expected_return, q{Added job_id to ALL } );

## Then job_id hash should stay the same except for ALL
# Clear ALL for this test
delete $job_id{ALL};

my @sample_job_ids = qw{ job_id_1 job_id_2 1234 };
is_deeply( \@{ $job_id{case1_MAIN}{sample1_MAIN} },
    \@sample_job_ids, q{Added dependencies to sample} );

done_testing();
