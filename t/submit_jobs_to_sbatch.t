#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Find;
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
use MIP::Constants qw{ $COMMA $COLON $EMPTY_STR $NEWLINE $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Processmanagement::Slurm_processes} => [qw{ submit_jobs_to_sbatch }],
        q{MIP::Test::Fixtures}                     => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Slurm_processes qw{ submit_jobs_to_sbatch };

diag(   q{Test submit_jobs_to_sbatch from Slurm_processes.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given a mock slurm instance
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );
my $sbatch_file_name =
  catfile( $Bin, qw{ data 643594-miptest test_script fastqc_ADM1059A1.0.sh } );

## Submit jobs to sbatch
my $job_id_returned = submit_jobs_to_sbatch(
    {
        base_command     => $slurm_mock_cmd,
        log              => $log,
        sbatch_file_name => $sbatch_file_name,
    }
);

## Then return a job id
my $expected_job_id = q{1234};
is( $job_id_returned, $expected_job_id, q{Returned expected job id} );

## Given a non integer return
my $slurm_mock_cmd_no_job_id =
  catfile( $Bin, qw{ data modules slurm-mock.pl } ) . q{ --no_job_id};

trap {
    submit_jobs_to_sbatch(
        {
            base_command     => $slurm_mock_cmd_no_job_id,
            log              => $log,
            sbatch_file_name => $sbatch_file_name,
        }
    )
};
## Required to make trap work - there must be something wrong with variables being stored in the trap object and not flushed between trappings
say {*STDOUT} $EMPTY_STR;

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if sbatch returned not a number} );
like(
    $trap->stderr,
    qr/Aborting \s+ run/xms,
    q{Throw fatal log message if sbatch returned not a number}
);

## Given a faulty sbatch script
trap {
    submit_jobs_to_sbatch(
        {
            log              => $log,
            sbatch_file_name => $sbatch_file_name,
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if faulty sbatch script} );

## Clean-up
# Needed to let slurm have time to fail
sleep 2 or croak( q{Cannot sleep } . $COLON . $OS_ERROR, $NEWLINE );

## Find file names in current working dir defined by code ref wanted
find( \&_wanted, cwd() );

sub _wanted {

    my $file_name = $_;    # Filename from find sub

    # It is not a file
    return if ( not -f $file_name );

    ## Only unlink slurm outfiles
    return if ( not $file_name =~ /slurm-\d+.out/sxm );

    unlink $file_name;
    return;
}

done_testing();
