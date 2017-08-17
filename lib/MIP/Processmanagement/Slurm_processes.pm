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
      qw(slurm_submit_job_no_dependency_add_to_case
	 slurm_submit_job_no_dependency_dead_end
	 slurm_submit_job_sample_id_dependency_dead_end
	 submit_jobs_to_sbatch slurm_submission_info);

}

## Constants
Readonly my $NEWLINE      => qq{\n};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $UNDERSCORE   => q{_};

sub slurm_submit_job_no_dependency_add_to_case {

##slurm_submit_job_no_dependency_add_to_case

##Function : Submit jobs that has no prior job dependencies and adds to case and samples dependencies using SLURM. Not dependent on earlier scripts and adds to sample_id jobs and family_id jobs, but sbatch is processed at family level i.e. affects all sample_id jobs e.g. building a reference (no_dependency_add_to_case).
##Returns  : ""
##Arguments: $job_id_href, $sample_ids_ref, $family_id, $path, $log, $sbatch_file_name
##         : $job_id_href      => Info on job id dependencies hash {REF}
##         : $sample_ids_ref   => Sample ids {REF}
##         : $family_id        => Family id
##         : $path             => Trunk or branch
##         : $log              => Log
##         : $sbatch_file_name => Sbatch file name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $sample_ids_ref;
    my $family_id;
    my $path;
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
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref
        },
        family_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
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

    use MIP::Processmanagement::Processes
      qw(add_sample_job_id_to_dependency_tree add_pan_job_id_to_family_id_dependency_tree);

    ## Set keys
    my $family_id_chain_key => $family_id . $UNDERSCORE . $path;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            sbatch_file_name => $sbatch_file_name,
            log              => $log,
        }
    );

  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        ## Add job_id_returned to hash

        # Set sample_id_chain_key
        my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

        ## Save job_id to sample_id dependencies
        add_sample_job_id_to_dependency_tree(
            {
                job_id_href         => $job_id_href,
                family_id_chain_key => $family_id_chain_key,
                sample_id_chain_key => $sample_id_chain_key,
                job_id_returned     => $job_id_returned,
            }
        );

        ## Save family_merged job_id to family dependencies
        add_pan_job_id_to_family_id_dependency_tree(
            {
                job_id_href         => $job_id_href,
                family_id_chain_key => $family_id_chain_key,
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

    ## Add job_id_returned to hash for sacct processing downstream
    push @{ $job_id_href->{PAN}{PAN} }, $job_id_returned;

    return;
}

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

    # The job_id that is returned from submission
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

sub slurm_submit_job_no_dependency_add_to_sample {

##slurm_submit_job_no_dependency_add_to_sample

##Function : Submit jobs that has no prior job dependencies but leave sample dependencies using SLURM
##Returns  : ""
##Arguments: $job_id_href, $family_id, $sample_id, $path, $sbatch_file_name, $log
##         : $job_id_href      => The info on job ids hash {REF}
##         : $family_id        => Family id
##         : $sample_id        => Sample id
##         : $path             => Trunk or branch
##         : $sbatch_file_name => Sbatch file name
##         : $log              => Log

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id;
    my $sample_id;
    my $path;
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
        family_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
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

    use MIP::Processmanagement::Processes
      qw(add_sample_job_id_to_dependency_tree);

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $family_id_chain_key => $family_id . $UNDERSCORE . $path;
    my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

    ## Initiate chain - No dependencies, initiate Trunk (Main or other)

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
            sbatch_file_name => $sbatch_file_name,
            log              => $log,
        }
    );

    add_sample_job_id_to_dependency_tree(
        {
            job_id_href         => $job_id_href,
            family_id_chain_key => $family_id_chain_key,
            sample_id_chain_key => $sample_id_chain_key,
            job_id_returned     => $job_id_returned,
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

sub slurm_submit_job_sample_id_dependency_dead_end {

##slurm_submit_job_sample_id_dependency_dead_end

##Function : Submit jobs that has sample_id dependencies and leave no dependencies using SLURM
##Returns  : ""
##Arguments: $job_id_href, $family_id, $sample_id, $path, $sbatch_file_name, $log, $job_dependency_type
##         : $job_id_href         => The info on job ids hash {REF}
##         : $family_id           => Family id
##         : $sample_id           => Sample id
##         : $path                => Trunk or branch
##         : $sbatch_file_name    => Sbatch file name
##         : $log                 => Log
##         : $job_dependency_type => SLURM job dependency type

    my ($arg_href) = @_;

    ## Default(s)
    my $job_dependency_type;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id;
    my $sample_id;
    my $path;
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
        family_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id
        },
        sample_id => {
            strict_type => 1,
            store       => \$sample_id
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
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
		job_dependency_type => {
            default     => q{afterok},
            allow       => [qw{afterany afterok}],
            strict_type => 1,
            store       => \$job_dependency_type
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    use MIP::Processmanagement::Processes
      qw(add_pan_job_id_to_sample_id_dependency_tree add_to_job_id_dependency_string);

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $family_id_chain_key = $family_id . $UNDERSCORE . $path;
    my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
    my $merged_chain_key = $family_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    ## Always check and add any pan (i.e job_ids that affect all chains) dependency jobs
      add_pan_job_id_to_sample_id_dependency_tree(
					   {
					    job_id_href         => $job_id_href,
					    family_id_chain_key => $family_id_chain_key,
					    sample_id_chain_key => $sample_id_chain_key,
					   }
					  );
      if ( $path eq q{MAIN} ) {

	## Add to job_id string
	$job_ids_string = add_to_job_id_dependency_string(
					{
					 job_id_href         => $job_id_href,
					 family_id_chain_key => $family_id_chain_key,
					 chain_key           => $sample_id_chain_key,
					}
				       );
      }
      elsif ( $path ne q{MAIN} ) {
	## Check for any previous job_ids within path current PATH. Branch.

	## Second or later in branch chain
	if ( $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} ) {

                ## Add to job_id string
                $job_ids_string = add_to_job_id_dependency_string(
                    {
                        job_id_href         => $job_id_href,
                        family_id_chain_key => $family_id_chain_key,
                        chain_key           => $sample_id_chain_key,
                    }
                );
            }
            elsif ( $job_id_href->{ $family_id . $UNDERSCORE . q{MAIN} }
                { $sample_id . $UNDERSCORE . q{MAIN} } )
            {
                ## Inherit from potential MAIN. Trunk

                ## Add to job_id string
                $job_ids_string = add_to_job_id_dependency_string(
                    {
                        job_id_href         => $job_id_href,
                        family_id_chain_key => $family_id . $UNDERSCORE . q{MAIN},
                        chain_key           => $sample_id . $UNDERSCORE . q{MAIN},
                    }
                );
            }
        }

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
	 job_dependency_type => $job_dependency_type,
	 job_ids_string      => $job_ids_string,
            sbatch_file_name => $sbatch_file_name,
            log              => $log,
        }
    );

            ## Clear latest family_id_sample_id chainkey
            @{ $job_id_href->{$family_id_chain_key}
                  { $merged_chain_key } } = ();

            #Clear latest sample_id chainkey
            @{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} } =
              ();

    ## Keeps the job_id string dependency within reasonable limits
        if (   ( defined $job_id_href->{ALL}{ALL} )
            && ( scalar( @{ $job_id_href->{ALL}{ALL} } ) >= 100 ) )
        {

            # Remove oldest job_id
            shift @{ $job_id_href->{ALL}{ALL} };
        }
        ## Job dependent on all jobs
        # Add job_id to hash
        push @{ $job_id_href->{ALL}{ALL} }, $job_id_returned;

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

        $log->fatal( @{$stderr_buf_ref} );
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


##Decommisioned###
sub slurm_submit_job_case_dependency_add_to_case {

##slurm_submit_job_case_dependency_add_to_case

##Function : Submit jobs that has case dependencies and leave case dependencies using SLURM
##Returns  : ""
##Arguments: $job_id_href, $family_id, $sample_id, $path, $sbatch_file_name, $log, $job_dependency_type
##         : $job_id_href         => The info on job ids hash {REF}
##         : $family_id           => Family id
##         : $sample_id           => Sample id
##         : $path                => Trunk or branch
##         : $sbatch_file_name    => Sbatch file name
##         : $log                 => Log
##         : $job_dependency_type => SLURM job dependency type

    my ($arg_href) = @_;

    ## Default(s)
    my $job_dependency_type;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id;
    my $sample_id;
    my $path;
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
        family_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id
        },
        sample_id => {
            strict_type => 1,
            store       => \$sample_id
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
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
		job_dependency_type => {
            default     => q{afterok},
            allow       => [qw{afterany afterok}],
            strict_type => 1,
            store       => \$job_dependency_type
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    use MIP::Processmanagement::Processes
      qw(add_pan_job_id_to_sample_id_dependency_tree);

    # Create string with all previous job_ids
    my $job_ids_string;

    # The job_id that is returned from submission
    my $job_id_returned;

    ## Set keys
    my $family_id_chain_key => $family_id . $UNDERSCORE . $path;

    if ( defined $sample_id ) {

      my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

      ## Saves job_id to the correct hash array depending on chaintype
      add_pan_job_id_to_sample_id_dependency_tree(
					   {
					    job_id_href         => $job_id_href,
					    family_id_chain_key => $family_id_chain_key,
					    sample_id_chain_key => $sample_id_chain_key,
					   }
					  );

      if ( $path eq q{MAIN} ) {

	## Add to job_id string
	$job_ids_string = add_to_job_id(
					{
					 job_id_href         => $job_id_href,
					 family_id_chain_key => $family_id_chain_key,
					 chain_key           => $sample_id_chain_key,
					}
				       );
      }
      elsif ( $path ne q{MAIN} ) {
	## Check for any previous job_ids within path current PATH. Branch.

	## Second or later in branch chain
	if ( $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} ) {

                ## Add to job_id string
                $job_ids_string = add_to_job_id(
                    {
                        job_id_href         => $job_id_href,
                        family_id_chain_key => $family_id_chain_key,
                        chain_key           => $sample_id_chain_key,
                    }
                );
            }
            elsif ( $job_id_href->{ $family_id . q{_MAIN} }
                { $sample_id . q{_MAIN} } )
            {
                ## Inherit from potential MAIN. Trunk

                ## Add to job_id string
                $job_ids_string = add_to_job_id(
                    {
                        job_id_href         => $job_id_href,
                        family_id_chain_key => $family_id . q{_MAIN},
                        chain_key           => $sample_id . q{_MAIN},
                    }
                );
            }
        }

    ## Submit jobs to sbatch
    $job_id_returned = submit_jobs_to_sbatch(
        {
	 job_dependency_type => $job_dependency_type,
	 job_ids_string      => $job_ids_string,
            sbatch_file_name => $sbatch_file_name,
            log              => $log,
        }
    );

            ## Clear latest family_id_sample_id chainkey
            @{ $job_id_href->{$family_id_chain_key}
                  { $family_id_chain_key . q{_} . $sample_id_chain_key } } = ();

            #Clear latest sample_id chainkey
            @{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} } =
              ();

            ## Add job_id_returned to hash
            push @{ $job_id_href->{$family_id_chain_key}
                  { $family_id_chain_key . q{_} . $sample_id_chain_key } },
              $job_id_returned;

    ## Keeps the job_id string dependecy within reasonable limits
        if (   ( defined $job_id_href->{ALL}{ALL} )
            && ( scalar( @{ $job_id_href->{ALL}{ALL} } ) >= 100 ) )
        {

            # Remove oldest job_id
            shift @{ $job_id_href->{ALL}{ALL} };
        }
        ## Job dependent on all jobs
        # Add job_id to hash
        push @{ $job_id_href->{ALL}{ALL} }, $job_id_returned;
    }
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

1;
