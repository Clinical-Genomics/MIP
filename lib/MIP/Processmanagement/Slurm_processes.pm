package MIP::Processmanagement::Slurm_processes;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    # Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];

use FindBin qw($Bin);    # Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );

BEGIN {
    use base qw (Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw(slurm_submit_job_no_dependency_dead_end submit_jobs_to_sbatch slurm_submission_info);

}

## Constants
Readonly my $NEWLINE      => qq{\n};
Readonly my $SINGLE_QUOTE => q{'};

sub slurm_submit_job_no_dependency_dead_end {

##slurm_submit_job_no_dependency_dead_end

##Function : Submit jobs that has no prior job dependencies and does not leave any dependencies using SLURM
##Returns  : ""
##Arguments: $job_id_href, $sbatch_file_name, $log
##         : $job_id_href      => The info on job ids hash {REF}
##         : $sbatch_file_name => Sbatch file name
##         : $log              => Log

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $sbatch_file_name;
    my $log;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        sbatch_file_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sbatch_file_name
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $job_id_returned;

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            sbatch_file_name => $sbatch_file_name,
            log              => $log,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add job_id_returned to hash for sacct processing downstream
    push @{ $job_id_href->{PAN}{PAN} }, $job_id_returned;

    return;
}

sub submit_jobs_to_sbatch {

##submit_sbatch_to_sbatch

##Function : Sumit jobs to sbatch
##Returns  : "The submitted $job_id"
##Arguments: $sbatch_file_name, $job_dependency_type, $job_ids_string, $log
##         : $sbatch_file_name    => Sbatch file to submit
##         : $job_dependency_type => Job dependency type
##         : $job_ids_string      => Job ids string
##         : $log                 => Log

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sbatch_file_name;
    my $log;
    my $job_dependency_type;
    my $job_ids_string;

    my $tmpl = {
        sbatch_file_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sbatch_file_name
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        job_dependency_type =>
          { strict_type => 1, store => \$job_dependency_type },
        job_ids_string => { strict_type => 1, store => \$job_ids_string },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $job_id;

    use MIP::Workloadmanager::Slurm qw(slurm_sbatch);
    use IPC::Cmd qw[ run ];

    # Supply with potential dependency of previous jobs that this one is dependent on
    my @commands = slurm_sbatch(
        {
            infile_path     => $sbatch_file_name,
            dependency_type => $job_dependency_type,
            job_ids_string  => $job_ids_string,
        }
    );

    # Submit job process
    my (
        $success_ref,    $error_message_ref, $full_buf_ref,
        $stdout_buf_ref, $stderr_buf_ref
      )
      = run(
        command => \@commands,
        verbose => 0
      );

    # Just submitted job_id
    if ( $stdout_buf_ref->[0] =~ /Submitted batch job (\d+)/ ) {

        $job_id = $1;
    }
    else {
       # Catch errors since, proper sbatch submission should only return numbers

        $log->fatal( @{ $stderr_buf_ref } );
        $log->fatal( q{Aborting run} . $NEWLINE );
        exit 1;
    }
    return $job_id;
}

sub slurm_submission_info {

##slurm_submission_info

##Function : Broadcast slurm submission info
##Returns  : ""
##Arguments: $log, $job_id_returned
##         : $log             => Log
##         : $job_id_returned => Job_id that was returned from submission

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $job_id_returned;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        job_id_returned => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$job_id_returned
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    $log->info( q{Sbatch script submitted, job id: } . $job_id_returned,
        $NEWLINE );
    $log->info(
        q{To check status of job, please run }
          . $SINGLE_QUOTE
          . q{squeue -j }
          . $job_id_returned
          . $SINGLE_QUOTE,
        $NEWLINE
    );
    $log->info(
        q{To cancel job, please run }
          . $SINGLE_QUOTE
          . q{scancel }
          . $job_id_returned
          . $SINGLE_QUOTE,
        $NEWLINE
    );
    return;
}

1;
