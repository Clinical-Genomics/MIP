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
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_log test_standard_cli };

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
        q{MIP::Processmanagement::Slurm_processes} =>
          [qw{ slurm_submit_job_no_dependency_dead_end }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Slurm_processes qw{ slurm_submit_job_no_dependency_dead_end };

diag(   q{Test slurm_submit_job_no_dependency_dead_end from Slurm_processes.pm v}
      . $MIP::Processmanagement::Slurm_processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a mock slurm and script
my %job_id         = test_mip_hashes( { mip_hash_name => q{job_id}, } );
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );
my $sbatch_file_name =
  catfile( $Bin, qw{ data 643594-miptest test_script fastqc_ADM1059A1.0.sh } );
my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

slurm_submit_job_no_dependency_dead_end(
    {
        base_command     => $slurm_mock_cmd,
        job_id_href      => \%job_id,
        log              => $log,
        sbatch_file_name => $sbatch_file_name,
    }
);

## Then add job_id returned to ALL
my $expected_return = q{1234};
is( $job_id{ALL}{ALL}[0], $expected_return, q{Added job_id to ALL } );

## Then job_id hash should stay the same except for pan
# Clear ALL for this test
delete $job_id{ALL};

my %original_job_id = test_mip_hashes( { mip_hash_name => q{job_id}, } );
is_deeply( \%job_id, \%original_job_id, q{Did not add to dependencies} );

done_testing();
