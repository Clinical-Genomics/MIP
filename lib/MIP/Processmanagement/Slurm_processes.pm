package MIP::Processmanagement::Slurm_processes;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };

BEGIN {
    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      slurm_submission_info
      slurm_submit_job_case_id_dependency_add_to_sample
      slurm_submit_chain_job_ids_dependency_add_to_path
      slurm_submit_job_no_dependency_dead_end
      slurm_submit_job_no_dependency_add_to_sample
      slurm_submit_job_no_dependency_add_to_samples
      slurm_submit_job_sample_id_dependency_add_to_sample
      slurm_submit_job_sample_id_dependency_add_to_case
      slurm_submit_job_sample_id_dependency_dead_end
      slurm_submit_job_sample_id_dependency_case_dead_end
      slurm_submit_job_sample_id_dependency_step_in_parallel
      slurm_submit_job_sample_id_dependency_step_in_parallel_to_case
      submit_jobs_to_sbatch
      submit_slurm_recipe
    };

}

sub slurm_submit_job_case_id_dependency_add_to_sample {

## Function : Submit jobs that has case_id dependencies and adds to sample dependencies using SLURM
## Returns  :
## Arguments: $base_command            => Sbatch
##          : $case_id                 => Case id
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_dependency_type     => SLURM job dependency type
##          : $job_id_href             => The info on job ids hash {REF}
##          : $log                     => Log
##          : $path                    => Trunk or branch
##          : $sample_id               => Sample id
##          : $sbatch_file_name        => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $path;
    my $sample_id;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;
    my $job_dependency_type;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_dependency_type => {
            allow       => [qw{ afterany afterok }],
            default     => q{afterok},
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        path      => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{
      add_job_id_dependency_tree
      add_pan_job_id_to_sample_id_dependency_tree
      add_parallel_job_id_to_sample_id_dependency_tree
      add_sample_job_id_to_sample_id_dependency_tree
      add_to_job_id_dependency_string
      create_job_id_string_for_case_id
      limit_job_id_string
    };

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $case_id_chain_key   = $case_id . $UNDERSCORE . $path;
    my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
    my $pan_chain_key       = $case_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    ## Always check and add any pan (i.e job_ids that affect all chains) dependency jobs
    add_pan_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    ## Always check and add any parallel (i.e job_ids that are processed in parallel witin path) dependency jobs
    add_parallel_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key       => $case_id_chain_key,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            path                    => $path,
            sample_id               => $sample_id,
            sample_id_chain_key     => $sample_id_chain_key,
        }
    );

    ## Create job id string from the job id chain and path associated with sample for
    ## SLURM submission using dependencies
    $job_ids_string = create_job_id_string_for_case_id(
        {
            case_id                 => $case_id,
            case_id_chain_key       => $case_id_chain_key,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            path                    => $path,
            sample_ids_ref          => [$sample_id],
        }
    );

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command        => $base_command,
            job_dependency_type => $job_dependency_type,
            job_ids_string      => $job_ids_string,
            log                 => $log,
            sbatch_file_name    => $sbatch_file_name,
        }
    );

    add_sample_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            job_id_returned     => $job_id_returned,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    ## Limit number of job_ids in job_id chain
    limit_job_id_string(
        {
            job_id_href => $job_id_href,
        }
    );

    ## Add job_id to jobs dependent on all jobs
    add_job_id_dependency_tree(
        {
            chain_key       => q{ALL},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_job_no_dependency_dead_end {

## Function : Submit jobs that has no prior job dependencies and does not leave any dependencies using SLURM, except to PAN dependencies
## Returns  :
## Arguments: $base_command     => Sbatch
##          : $log              => Log object
##          : $job_id_href      => The info on job ids hash {REF}
##          : $sbatch_file_name => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $job_id_href;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{ add_job_id_dependency_tree };

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command     => $base_command,
            log              => $log,
            sbatch_file_name => $sbatch_file_name,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_job_no_dependency_add_to_sample {

## Function : Submit jobs that has no prior job dependencies but leave sample dependencies using SLURM
## Returns  :
## Arguments: $base_command     => Sbatch
##          : $case_id          => Case id
##          : $job_id_href      => The info on job ids hash {REF}
##          : $log              => Log
##          : $path             => Trunk or branch
##          : $sample_id        => Sample id
##          : $sbatch_file_name => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $job_id_href;
    my $log;
    my $path;
    my $sample_id;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        path      => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes
      qw{add_sample_job_id_to_sample_id_dependency_tree
      add_job_id_dependency_tree};

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $case_id_chain_key   = $case_id . $UNDERSCORE . $path;
    my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

    ## Initiate chain - No dependencies, initiate Trunk (Main or other)

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command     => $base_command,
            log              => $log,
            sbatch_file_name => $sbatch_file_name,
        }
    );

    add_sample_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            job_id_returned     => $job_id_returned,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_job_no_dependency_add_to_samples {

## Function : Submit jobs that has no prior job dependencies and adds to samples dependencies using SLURM. Not dependent on earlier scripts and adds to sample_id jobs, but sbatch is processed at case level i.e. affects all sample_id jobs e.g. building a reference.
## Returns  :
## Arguments: $base_command     => Sbatch
##          : $case_id          => Case id
##          : $job_id_href      => Info on job id dependencies hash {REF}
##          : $log              => Log
##          : $path             => Trunk or branch
##          : $sample_ids_ref   => Sample ids {REF}
##          : $sbatch_file_name => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $job_id_href;
    my $log;
    my $path;
    my $sample_ids_ref;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes
      qw{add_sample_job_id_to_sample_id_dependency_tree add_job_id_dependency_tree};

    ## Set keys
    my $case_id_chain_key = $case_id . $UNDERSCORE . $path;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command     => $base_command,
            log              => $log,
            sbatch_file_name => $sbatch_file_name,
        }
    );

  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        ## Add job_id_returned to hash

        # Set sample_id_chain_key
        my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

        ## Save job_id to sample_id dependencies
        add_sample_job_id_to_sample_id_dependency_tree(
            {
                case_id_chain_key   => $case_id_chain_key,
                job_id_href         => $job_id_href,
                job_id_returned     => $job_id_returned,
                sample_id_chain_key => $sample_id_chain_key,
            }
        );
    }

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_job_sample_id_dependency_dead_end {

## Function : Submit jobs that has sample_id dependencies and leave no dependencies using SLURM
## Returns  :
## Arguments: $base_command            => Sbatch
##          : $case_id                 => Case id
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => The info on job ids hash {REF}
##          : $job_dependency_type     => SLURM job dependency type
##          : $log                     => Log
##          : $path                    => Trunk or branch
##          : $sample_id               => Sample id
##          : $sbatch_file_name        => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $path;
    my $sample_id;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;
    my $job_dependency_type;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_dependency_type => {
            allow       => [qw{afterany afterok}],
            default     => q{afterok},
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        path      => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{add_job_id_dependency_tree
      add_pan_job_id_to_sample_id_dependency_tree
      add_parallel_job_id_to_sample_id_dependency_tree
      add_to_job_id_dependency_string
      create_job_id_string_for_sample_id
      clear_pan_job_id_dependency_tree
      clear_sample_id_job_id_dependency_tree
      clear_sample_id_parallel_job_id_dependency_tree
      limit_job_id_string
    };

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $case_id_chain_key   = $case_id . $UNDERSCORE . $path;
    my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
    my $pan_chain_key       = $case_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    ## Always check and add any pan (i.e job_ids that affect all chains) dependency jobs
    add_pan_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    ## Always check and add any parallel (i.e job_ids that are processed in parallel witin path) dependency jobs
    add_parallel_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key       => $case_id_chain_key,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            path                    => $path,
            sample_id               => $sample_id,
            sample_id_chain_key     => $sample_id_chain_key,
        }
    );

    ## Create job id string from the job id chain and path associated with sample for
    ## SLURM submission using dependencies
    $job_ids_string = create_job_id_string_for_sample_id(
        {
            case_id             => $case_id,
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            path                => $path,
            sample_id           => $sample_id,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command        => $base_command,
            job_dependency_type => $job_dependency_type,
            job_ids_string      => $job_ids_string,
            log                 => $log,
            sbatch_file_name    => $sbatch_file_name,
        }
    );

    ## Clear latest case_id_sample_id chainkey
    clear_pan_job_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    clear_sample_id_parallel_job_id_dependency_tree(
        {
            case_id_chain_key       => $case_id_chain_key,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            path                    => $path,
            sample_id               => $sample_id,
        }
    );
    ## Clear sample job_ids in the the sample_id chain
    clear_sample_id_job_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    ## Limit number of job_ids in job_id chain
    limit_job_id_string(
        {
            job_id_href => $job_id_href,
        }
    );

    ## Add job_id to jobs dependent on all jobs
    add_job_id_dependency_tree(
        {
            chain_key       => q{ALL},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_job_sample_id_dependency_add_to_sample {

## Function : Submit jobs that has sample_id dependencies and adds to sample dependencies using SLURM
## Returns  :
## Arguments: $base_command            => Sbatch
##          : $case_id                 => Case id
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_dependency_type     => SLURM job dependency type
##          : $job_id_href             => The info on job ids hash {REF}
##          : $log                     => Log
##          : $path                    => Trunk or branch
##          : $sample_id               => Sample id
##          : $sbatch_file_name        => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $path;
    my $sample_id;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;
    my $job_dependency_type;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_dependency_type => {
            allow       => [qw{ afterany afterok }],
            default     => q{afterok},
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        path      => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{
      add_job_id_dependency_tree
      add_pan_job_id_to_sample_id_dependency_tree
      add_to_job_id_dependency_string
      clear_pan_job_id_dependency_tree
      clear_sample_id_job_id_dependency_tree
      clear_sample_id_parallel_job_id_dependency_tree
      create_job_id_string_for_sample_id
      limit_job_id_string
    };

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $case_id_chain_key   = $case_id . $UNDERSCORE . $path;
    my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
    my $pan_chain_key       = $case_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    ## Always check and add any pan (i.e job_ids that affect all chains) dependency jobs
    add_pan_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    ## Always check and add any parallel (i.e job_ids that are processed in parallel witin path) dependency jobs
    add_parallel_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key       => $case_id_chain_key,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            path                    => $path,
            sample_id               => $sample_id,
            sample_id_chain_key     => $sample_id_chain_key,
        }
    );

    ## Create job id string from the job id chain and path associated with sample for
    ## SLURM submission using dependencies
    $job_ids_string = create_job_id_string_for_sample_id(
        {
            case_id             => $case_id,
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            path                => $path,
            sample_id           => $sample_id,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command        => $base_command,
            job_dependency_type => $job_dependency_type,
            job_ids_string      => $job_ids_string,
            log                 => $log,
            sbatch_file_name    => $sbatch_file_name,
        }
    );

    ## Clear latest case_id_sample_id chainkey
    clear_pan_job_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    clear_sample_id_parallel_job_id_dependency_tree(
        {
            case_id_chain_key       => $case_id_chain_key,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            path                    => $path,
            sample_id               => $sample_id,
        }
    );

    ## Clear sample job_ids in the the sample_id chain
    clear_sample_id_job_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    add_sample_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            job_id_returned     => $job_id_returned,
            sample_id_chain_key => $sample_id_chain_key,
        }
    );

    ## Limit number of job_ids in job_id chain
    limit_job_id_string(
        {
            job_id_href => $job_id_href,
        }
    );

    ## Add job_id to jobs dependent on all jobs
    add_job_id_dependency_tree(
        {
            chain_key       => q{ALL},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_job_sample_id_dependency_add_to_case {

## Function : Submit jobs that has sample_id dependencies and adds to case dependencies using SLURM
## Returns  :
## Arguments: $base_command            => Sbatch
##          : $case_id                 => Case id
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_dependency_type     => SLURM job dependency type
##          : $job_id_href             => The info on job ids hash {REF}
##          : $log                     => Log
##          : $parallel_chains_ref     => Info on parallel chains array {REF}
##          : $path                    => Trunk or branch
##          : $sample_ids_ref          => Sample ids {REF}
##          : $sbatch_file_name        => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $parallel_chains_ref;
    my $path;
    my $sample_ids_ref;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;
    my $job_dependency_type;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_dependency_type => {
            allow       => [qw{ afterany afterok }],
            default     => q{afterok},
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parallel_chains_ref => {
            default     => [],
            store       => \$parallel_chains_ref,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{
      add_job_id_dependency_tree
      add_parallel_job_ids_to_job_id_dependency_string
      clear_all_job_ids_within_chain_key_dependency_tree
      clear_case_id_job_id_dependency_tree
      create_job_id_string_for_case_id
      limit_job_id_string
    };

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $case_id_chain_key = $case_id . $UNDERSCORE . $path;

    ## Create job id string from the job id chain and path associated with sample for
    ## SLURM submission using dependencies
    $job_ids_string = create_job_id_string_for_case_id(
        {
            case_id                 => $case_id,
            case_id_chain_key       => $case_id_chain_key,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            parallel_chains_ref     => $parallel_chains_ref,
            path                    => $path,
            sample_ids_ref          => $sample_ids_ref,
        }
    );

    ## Check if last step submission was parallel
    my $parallel_job_ids_string = add_parallel_job_ids_to_job_id_dependency_string(
        {
            case_id_chain_key => $case_id_chain_key,
            job_id_href       => $job_id_href,
        }
    );

    ## If parellel job_ids existed add to current job_id_string
    if ($parallel_job_ids_string) {

        $job_ids_string .= $parallel_job_ids_string;
    }

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command        => $base_command,
            job_dependency_type => $job_dependency_type,
            job_ids_string      => $job_ids_string,
            log                 => $log,
            sbatch_file_name    => $sbatch_file_name,
        }
    );

    ## Clear all job_ids in the the chain key
    clear_all_job_ids_within_chain_key_dependency_tree(
        {
            case_id_chain_key => $case_id_chain_key,
            job_id_href       => $job_id_href,
        }
    );

    ## Add case_id job_id to case dependency tree
    add_job_id_dependency_tree(
        {
            chain_key       => $case_id_chain_key,
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );

    ## Limit number of job_ids in job_id chain
    limit_job_id_string(
        {
            job_id_href => $job_id_href,
        }
    );

    ## Add job_id to jobs dependent on all jobs
    add_job_id_dependency_tree(
        {
            chain_key       => q{ALL},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_job_sample_id_dependency_case_dead_end {

## Function : Submit jobs that has sample_id dependencies and leave no dependencies using SLURM
## Returns  :
## Arguments: $base_command            => Sbatch
##          : $case_id                 => Case id
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_dependency_type     => SLURM job dependency type
##          : $job_id_href             => The info on job ids hash {REF}
##          : $log                     => Log
##          : $parallel_chains_ref     => Info on parallel chains array {REF}
##          : $path                    => Trunk or branch
##          : $sample_ids_ref          => Sample ids {REF}
##          : $sbatch_file_name        => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $parallel_chains_ref;
    my $path;
    my $sample_ids_ref;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;
    my $job_dependency_type;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        job_dependency_type => {
            allow       => [qw{ afterany afterok }],
            default     => q{afterok},
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parallel_chains_ref => {
            default     => [],
            store       => \$parallel_chains_ref,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{
      add_job_id_dependency_tree
      add_parallel_job_ids_to_job_id_dependency_string
      clear_all_job_ids_within_chain_key_dependency_tree
      clear_case_id_job_id_dependency_tree
      create_job_id_string_for_case_id
      limit_job_id_string
    };

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $case_id_chain_key = $case_id . $UNDERSCORE . $path;

    ## Create job id string from the job id chain and path associated with sample for
    ## SLURM submission using dependencies
    $job_ids_string = create_job_id_string_for_case_id(
        {
            case_id                 => $case_id,
            case_id_chain_key       => $case_id_chain_key,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            parallel_chains_ref     => $parallel_chains_ref,
            path                    => $path,
            sample_ids_ref          => $sample_ids_ref,
        }
    );

    ## Check if last case job submission was parallel
    my $parallel_job_ids_string = add_parallel_job_ids_to_job_id_dependency_string(
        {
            case_id_chain_key => $case_id_chain_key,
            job_id_href       => $job_id_href,
        }
    );

    ## If parellel job_ids existed add to current job_id_string
    if ($parallel_job_ids_string) {

        $job_ids_string .= $parallel_job_ids_string;
    }

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command        => $base_command,
            job_dependency_type => $job_dependency_type,
            job_ids_string      => $job_ids_string,
            log                 => $log,
            sbatch_file_name    => $sbatch_file_name,
        }
    );

    ## Clear all job_ids in the the chain key
    clear_all_job_ids_within_chain_key_dependency_tree(
        {
            case_id_chain_key => $case_id_chain_key,
            job_id_href       => $job_id_href,
        }
    );

    ## Limit number of job_ids in job_id chain
    limit_job_id_string(
        {
            job_id_href => $job_id_href,
        }
    );

    ## Add job_id to jobs dependent on all jobs
    add_job_id_dependency_tree(
        {
            chain_key       => q{ALL},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_job_sample_id_dependency_step_in_parallel_to_case {

## Function : Submit jobs that has sample_id dependencies and adds to case parallel dependencies using SLURM
## Returns  :
## Arguments: $base_command            => Sbatch
##          : $case_id                 => Case id
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_dependency_type     => SLURM job dependency type
##          : $job_id_href             => The info on job ids hash {REF}
##          : $log                     => Log
##          : $parallel_chains_ref     => Info on parallel chains array {REF}
##          : $path                    => Trunk or branch
##          : $sample_ids_ref          => Sample ids {REF}
##          : $sbatch_file_name        => Sbatch file name
##          : $sbatch_script_tracker   => Track the number of parallel processes (e.g. sbatch scripts for a module)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $parallel_chains_ref;
    my $path;
    my $sample_ids_ref;
    my $sbatch_file_name;
    my $sbatch_script_tracker;

    ## Default(s)
    my $base_command;
    my $job_dependency_type;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_dependency_type => {
            allow       => [qw{afterany afterok}],
            default     => q{afterok},
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        parallel_chains_ref => {
            default     => [],
            store       => \$parallel_chains_ref,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
        sbatch_script_tracker => {
            allow       => qr{ \A\d+\z }sxm,
            defined     => 1,
            required    => 1,
            store       => \$sbatch_script_tracker,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{
      add_job_id_dependency_tree
      add_parallel_job_id_to_parallel_dependency_tree
      clear_case_id_job_id_dependency_tree
      create_job_id_string_for_case_id
      limit_job_id_string
    };

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $case_id_chain_key = $case_id . $UNDERSCORE . $path;

    ## Create job id string from the job id chain and path associated with sample for
    ## SLURM submission using dependencies
    $job_ids_string = create_job_id_string_for_case_id(
        {
            case_id                 => $case_id,
            case_id_chain_key       => $case_id_chain_key,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            parallel_chains_ref     => $parallel_chains_ref,
            path                    => $path,
            sample_ids_ref          => $sample_ids_ref,
        }
    );

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command        => $base_command,
            job_dependency_type => $job_dependency_type,
            job_ids_string      => $job_ids_string,
            log                 => $log,
            sbatch_file_name    => $sbatch_file_name,
        }
    );

    ## Add parallel case_id job_id to case dependency tree
    add_parallel_job_id_to_parallel_dependency_tree(
        {
            case_id_chain_key     => $case_id_chain_key,
            id                    => $case_id,
            job_id_href           => $job_id_href,
            job_id_returned       => $job_id_returned,
            path                  => $path,
            sbatch_script_tracker => $sbatch_script_tracker,
        }
    );

    ## Limit number of job_ids in job_id chain
    limit_job_id_string(
        {
            job_id_href => $job_id_href,
        }
    );

    ## Add job_id to jobs dependent on all jobs
    add_job_id_dependency_tree(
        {
            chain_key       => q{ALL},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_job_sample_id_dependency_step_in_parallel {

## Function : Submit jobs that has sample_id dependencies and are processed in parallel dependencies using SLURM
## Returns  :
## Arguments: $base_command            => Sbatch
##          : $case_id                 => Case id
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_dependency_type     => SLURM job dependency type
##          : $job_id_href             => The info on job ids hash {REF}
##          : $log                     => Log
##          : $path                    => Trunk or branch
##          : $sample_id               => Sample id
##          : $sbatch_file_name        => Sbatch file name
##          : $sbatch_script_tracker   => Track the number of parallel processes (e.g. sbatch scripts for a module)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $path;
    my $sample_id;
    my $sbatch_file_name;
    my $sbatch_script_tracker;

    ## Default(s)
    my $base_command;
    my $job_dependency_type;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        job_dependency_type => {
            allow       => [qw{ afterany afterok }],
            default     => q{afterok},
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        path      => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
        sbatch_script_tracker => {
            allow       => qr{ \A\d+\z }sxm,
            defined     => 1,
            required    => 1,
            store       => \$sbatch_script_tracker,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{
      add_job_id_dependency_tree
      add_pan_job_id_to_sample_id_dependency_tree
      add_sample_job_id_to_sample_id_dependency_tree
      add_to_job_id_dependency_string
      clear_pan_job_id_dependency_tree
      clear_sample_id_job_id_dependency_tree
      clear_sample_id_parallel_job_id_dependency_tree
      create_job_id_string_for_sample_id
      limit_job_id_string
    };

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $case_id_chain_key   = $case_id . $UNDERSCORE . $path;
    my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
    my $pan_chain_key       = $case_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    # Sample parallel chainkey
    my $sample_id_parallel_chain_key =
        $sample_id
      . $UNDERSCORE
      . q{parallel}
      . $UNDERSCORE
      . $path
      . $sbatch_script_tracker;

    ### For parallel steps the sub will be used multiple times so all prior job_id
    ### dependecies has to be pushed to the job_id_string directly and not added to
    ### the sample_id dependency tree

    ## Add prior parallel job_ids for this sbatch script tracker number to job_id string
    my $parallel_job_ids_string = add_to_job_id_dependency_string(
        {
            chain_key         => $sample_id_parallel_chain_key,
            case_id_chain_key => $case_id_chain_key,
            job_id_href       => $job_id_href,
        }
    );

    ## Add prior pan job_ids to job_id string
    my $pan_job_ids_string = add_to_job_id_dependency_string(
        {
            chain_key         => $pan_chain_key,
            case_id_chain_key => $case_id_chain_key,
            job_id_href       => $job_id_href,
        }
    );

    ## Create job id string from the job id chain and path associated with sample for
    ## SLURM submission using dependencies
    $job_ids_string = create_job_id_string_for_sample_id(
        {
            case_id               => $case_id,
            case_id_chain_key     => $case_id_chain_key,
            job_id_href           => $job_id_href,
            path                  => $path,
            sample_id             => $sample_id,
            sample_id_chain_key   => $sample_id_chain_key,
            sbatch_script_tracker => $sbatch_script_tracker,
        }
    );

    ## Concatenate other job id strings to current job id string
    my @job_id_strings = ( $parallel_job_ids_string, $pan_job_ids_string );

    foreach my $other_job_id_string (@job_id_strings) {

        if ($other_job_id_string) {

            $job_ids_string .= $other_job_id_string;
        }
    }

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command        => $base_command,
            job_dependency_type => $job_dependency_type,
            job_ids_string      => $job_ids_string,
            log                 => $log,
            sbatch_file_name    => $sbatch_file_name,
        }
    );

    ## Add to parallel chain dependency for this sbatch script tracker number
    add_sample_job_id_to_sample_id_dependency_tree(
        {
            case_id_chain_key   => $case_id_chain_key,
            job_id_href         => $job_id_href,
            job_id_returned     => $job_id_returned,
            sample_id_chain_key => $sample_id_parallel_chain_key,
        }
    );

    ## Limit number of job_ids in job_id chain
    limit_job_id_string(
        {
            job_id_href => $job_id_href,
        }
    );

    ## Add job_id to jobs dependent on all jobs
    add_job_id_dependency_tree(
        {
            chain_key       => q{ALL},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub slurm_submit_chain_job_ids_dependency_add_to_path {

## Function : Submit jobs that has all dependencies and adds to all dependencies using SLURM
## Returns  :
## Arguments: $base_command        => Sbatch
##          : $job_dependency_type => SLURM job dependency type
##          : $job_id_href         => The info on job ids hash {REF}
##          : $log                 => Log
##          : $path                => Trunk or branch
##          : $sbatch_file_name    => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $log;
    my $path;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;
    my $job_dependency_type;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        job_dependency_type => {
            allow       => [qw{ afterany afterok }],
            default     => q{afterok},
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{
      add_job_id_dependency_tree
      add_to_job_id_dependency_string
      limit_job_id_string
    };

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $case_id_chain_key = $path;
    my $chain_key         = $path;

    $job_ids_string = add_to_job_id_dependency_string(
        {
            chain_key         => $chain_key,
            case_id_chain_key => $case_id_chain_key,
            job_id_href       => $job_id_href,
        }
    );

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            base_command        => $base_command,
            job_dependency_type => $job_dependency_type,
            job_ids_string      => $job_ids_string,
            log                 => $log,
            sbatch_file_name    => $sbatch_file_name,
        }
    );

    ## Limit number of job_ids in job_id chain
    limit_job_id_string(
        {
            job_id_href => $job_id_href,
        }
    );

    ## Add job_id to jobs dependent on all jobs
    add_job_id_dependency_tree(
        {
            chain_key       => $chain_key,
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );

    slurm_submission_info(
        {
            job_id_returned => $job_id_returned,
            log             => $log,
        }
    );

    ## Add PAN job_id_returned to hash for sacct processing downstream
    add_job_id_dependency_tree(
        {
            chain_key       => q{PAN},
            job_id_href     => $job_id_href,
            job_id_returned => $job_id_returned,
        }
    );
    return;
}

sub submit_jobs_to_sbatch {

## Function : Sumit jobs to sbatch using dependencies
## Returns  : $job_id
## Arguments: $base_command        => Sbatch
##          : $job_dependency_type => Job dependency type
##          : $job_ids_string      => Job ids string
##          : $log                 => Log
##          : $sbatch_file_name    => Sbatch file to submit

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_dependency_type;
    my $job_ids_string;
    my $log;
    my $sbatch_file_name;

    ## Default(s)
    my $base_command;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        job_dependency_type => { store => \$job_dependency_type, strict_type => 1, },
        job_ids_string      => { store => \$job_ids_string,      strict_type => 1, },
        log                 => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        sbatch_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$sbatch_file_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Workloadmanager::Slurm qw{slurm_sbatch};
    use MIP::Unix::System qw{ system_cmd_call };

    ## Supply with potential dependency of previous jobs that this one is dependent on
    my @commands = slurm_sbatch(
        {
            base_command    => $base_command,
            dependency_type => $job_dependency_type,
            infile_path     => $sbatch_file_name,
            job_ids_string  => $job_ids_string,
        }
    );

    # Submit job process
    my %return = system_cmd_call( { command_string => join $SPACE, @commands, } );

    # Sbatch should return message and job id in stdout
    croak(
        $log->fatal( @{ $return{error} } ) . $log->fatal( q{Aborting run} . $NEWLINE ) )
      if ( not $return{output}[0] );

    # Capture job id for submitted scripts
    my ($job_id) = $return{output}[0] =~ /Submitted \s+ batch \s+ job \s+ (\d+)/sxm;

    # Sbatch should return message and job id in stdout
    croak(
        $log->fatal( @{ $return{error} } ) . $log->fatal( q{Aborting run} . $NEWLINE ) )
      if ( not $job_id );

    return $job_id;
}

sub slurm_submission_info {

## Function : Broadcast slurm submission info
## Returns  :
## Arguments: $job_id_returned => Job_id that was returned from submission
##          : $log             => Log

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_returned;
    my $log;

    my $tmpl = {
        job_id_returned => {
            defined     => 1,
            required    => 1,
            store       => \$job_id_returned,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $log->info( q{Sbatch script submitted, job id: } . $job_id_returned, $NEWLINE );
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

sub submit_slurm_recipe {

## Function : Submit SLURM recipe
## Returns  :
## Arguments: $base_command            => Sbatch
##          : $case_id                 => Case id
##          : $dependency_method       => Dependency method
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_dependency_type     => SLURM job dependency type
##          : $job_id_chain            => Chain id
##          : $job_id_href             => The info on job ids hash {REF}
##          : $job_dependency_type     => SLURM job dependency type
##          : $log                     => Log
##          : $parallel_chains_ref     => Info on parallel chains array {REF}
##          : $sample_id               => Sample id
##          : $sample_ids_ref          => Sample ids {REF}
##          : $recipe_file_path        => Sbatch file path
##          : $recipe_files_tracker    => Track the number of parallel processes (e.g. sbatch scripts for a module)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dependency_method;
    my $case_id;
    my $infile_lane_prefix_href;
    my $job_dependency_type;
    my $job_id_chain;
    my $job_id_href;
    my $log;
    my $parallel_chains_ref;
    my $sample_id;
    my $sample_ids_ref;
    my $recipe_file_path;
    my $recipe_files_tracker;

    ## Default(s)
    my $base_command;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        case_id => {
            store       => \$case_id,
            strict_type => 1,
        },
        dependency_method => {
            allow => [
                qw{ add_to_all
                  case_to_island
                  island
                  island_to_sample
                  island_to_samples
                  case_to_sample
                  sample_to_case
                  sample_to_case_parallel
                  sample_to_island
                  sample_to_sample
                  sample_to_sample_parallel }
            ],
            defined     => 1,
            required    => 1,
            store       => \$dependency_method,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_dependency_type => {
            allow       => [ undef, qw{ afterany afterok } ],
            default     => q{afterok},
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        job_id_chain => { store => \$job_id_chain, strict_type => 1, },
        job_id_href  => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parallel_chains_ref => {
            default     => [],
            store       => \$parallel_chains_ref,
            strict_type => 1,
        },
        recipe_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_file_path,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        recipe_files_tracker => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            store       => \$recipe_files_tracker,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $dependency_method eq q{add_to_all} ) {

        slurm_submit_chain_job_ids_dependency_add_to_path(
            {
                base_command        => $base_command,
                job_id_href         => $job_id_href,
                job_dependency_type => $job_dependency_type,
                log                 => $log,
                path                => $job_id_chain,
                sbatch_file_name    => $recipe_file_path,
            }
        );

        return 1;
    }
    if ( $dependency_method eq q{case_to_island} ) {

        slurm_submit_job_sample_id_dependency_case_dead_end(
            {
                base_command            => $base_command,
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref          => $sample_ids_ref,
                sbatch_file_name        => $recipe_file_path,
            }
        );
        return 1;
    }
    if ( $dependency_method eq q{island} ) {

        ## No upstream or downstream dependencies
        slurm_submit_job_no_dependency_dead_end(
            {
                base_command     => $base_command,
                job_id_href      => $job_id_href,
                log              => $log,
                sbatch_file_name => $recipe_file_path,
            }
        );
        return 1;
    }
    if ( $dependency_method eq q{island_to_sample} ) {

        slurm_submit_job_no_dependency_add_to_sample(
            {
                base_command     => $base_command,
                case_id          => $case_id,
                job_id_href      => $job_id_href,
                log              => $log,
                path             => $job_id_chain,
                sample_id        => $sample_id,
                sbatch_file_name => $recipe_file_path
            }
        );
        return 1;
    }
    if ( $dependency_method eq q{island_to_samples} ) {

        slurm_submit_job_no_dependency_add_to_samples(
            {
                base_command     => $base_command,
                case_id          => $case_id,
                job_id_href      => $job_id_href,
                log              => $log,
                path             => $job_id_chain,
                sample_ids_ref   => $sample_ids_ref,
                sbatch_file_name => $recipe_file_path,
            }
        );
        return 1;
    }
    if ( $dependency_method eq q{case_to_sample} ) {

        slurm_submit_job_case_id_dependency_add_to_sample(
            {
                base_command            => $base_command,
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $recipe_file_path,
            }
        );

    }
    if ( $dependency_method eq q{sample_to_case} ) {

        slurm_submit_job_sample_id_dependency_add_to_case(
            {
                base_command            => $base_command,
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                parallel_chains_ref     => $parallel_chains_ref,
                sample_ids_ref          => $sample_ids_ref,
                sbatch_file_name        => $recipe_file_path,
            }
        );
        return 1;
    }
    if ( $dependency_method eq q{sample_to_case_parallel} ) {

        slurm_submit_job_sample_id_dependency_step_in_parallel_to_case(
            {
                base_command            => $base_command,
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref          => $sample_ids_ref,
                sbatch_file_name        => $recipe_file_path,
                sbatch_script_tracker   => $recipe_files_tracker,
            }
        );
        return 1;
    }
    if ( $dependency_method eq q{sample_to_island} ) {

        slurm_submit_job_sample_id_dependency_dead_end(
            {
                base_command            => $base_command,
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $recipe_file_path,
            }
        );
        return 1;
    }
    if ( $dependency_method eq q{sample_to_sample_parallel} ) {

        slurm_submit_job_sample_id_dependency_step_in_parallel(
            {
                base_command            => $base_command,
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $recipe_file_path,
                sbatch_script_tracker   => $recipe_files_tracker,
            }
        );
        return 1;
    }
    ## Last dependency method
    return if ( $dependency_method ne q{sample_to_sample} );

    slurm_submit_job_sample_id_dependency_add_to_sample(
        {
            base_command            => $base_command,
            case_id                 => $case_id,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            log                     => $log,
            path                    => $job_id_chain,
            sample_id               => $sample_id,
            sbatch_file_name        => $recipe_file_path,
        }
    );
    return 1;
}

1;
