package MIP::Processmanagement::Processes;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile };
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
use MIP::Constants qw{ $DOT $EMPTY_STR $COLON $LOG_NAME $NEWLINE $UNDERSCORE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      add_to_job_id_dependency_string
      add_sample_ids_job_ids_to_job_id_dependency_string
      add_parallel_job_ids_to_job_id_dependency_string
      add_parallel_chains_job_ids_to_job_id_dependency_string
      add_job_id_dependency_tree
      add_parallel_job_id_to_sample_id_dependency_tree
      add_parallel_job_id_to_case_id_dependency_tree
      add_parallel_job_id_to_parallel_dependency_tree
      add_sample_id_parallel_job_id_to_case_id_dependency_tree
      add_sample_ids_parallel_job_id_to_case_id_dependency_tree
      add_pan_job_id_to_sample_id_dependency_tree
      add_sample_job_id_to_sample_id_dependency_tree
      add_sample_job_id_to_case_id_dependency_tree
      create_job_id_string_for_sample_id
      create_job_id_string_for_case_id
      create_job_id_string_for_case_id_and_path
      clear_sample_id_parallel_job_id_dependency_tree
      clear_pan_job_id_dependency_tree
      clear_sample_id_job_id_dependency_tree
      clear_case_id_job_id_dependency_tree
      clear_all_job_ids_within_chain_key_dependency_tree
      get_all_job_ids
      limit_job_id_string
      print_wait
      submit_recipe
      write_job_ids_to_file
    };
}

sub add_to_job_id_dependency_string {

## Function : Adds all previous job_ids per case_chain_key and chain_key to job_ids dependency string, which is used to set the dependency in SLURM.
## Returns  : $job_ids
## Arguments: $chain_key         => The current chain hash key
##          : $job_id_href       => Info on jobIds hash {REF}
##          : $case_id_chain_key => Case id chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_key;
    my $case_id_chain_key;
    my $job_id_href;

    my $tmpl = {
        chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$chain_key,
            strict_type => 1,
        },
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Job_id string to submit to workload manager
    my $job_ids_string = $EMPTY_STR;

    # Alias hash to get job_id dependencies from
    my $chain_job_id_href = $job_id_href->{$case_id_chain_key}{$chain_key};

    if ($chain_job_id_href) {

      JOB_IDS:
        while ( my ( $job_index, $job_id ) = each @{$chain_job_id_href} ) {

            # Only for defined job_ids
            if ( defined $job_id ) {

                ## Only 1 previous job_id
                if ( ( !$job_index ) && ( scalar @{$chain_job_id_href} == 1 ) ) {

                    # Single job_id start with $COLON and end without $COLON
                    $job_ids_string .= $COLON . $job_id;
                }
                elsif ( !$job_index ) {
                    ## First job_id to add

                    # First job_id start with $COLON and ends with $COLON
                    $job_ids_string .= $COLON . $job_id . $COLON;
                }
                elsif ( $job_index eq ( scalar @{$chain_job_id_href} - 1 ) ) {
                    ## Last job_id

                    # Last job_id finish without $COLON
                    $job_ids_string .= $job_id;
                }
                else {
                    ## JobIDs in the middle

                    $job_ids_string .= $job_id . $COLON;
                }
            }
        }
    }
    return $job_ids_string;
}

sub add_sample_ids_job_ids_to_job_id_dependency_string {

## Function : Create job id string from sample_ids job id chain and path for SLURM submission using dependencies
## Returns  : $job_ids_string
## Arguments: $case_id                           => Case id
##          : $case_id_chain_key                 => Case id chain hash key
##          : $job_id_href                       => The info on job ids hash {REF}
##          : $max_parallel_processes_count_href => Maximum number of parallel processes
##          : $path                              => Trunk or branch
##          : $sample_ids_ref                    => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $case_id_chain_key;
    my $job_id_href;
    my $max_parallel_processes_count_href;
    my $path;
    my $sample_ids_ref;

    my $tmpl = {
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        max_parallel_processes_count_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$max_parallel_processes_count_href,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;

    ## Add all previous jobId(s) from sample_id chain_key(s)
  SAMPLE_IDS:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

        if ( $job_id_href->{$case_id_chain_key}{$sample_id_chain_key} ) {

            ## Add to job_id string
            $job_ids_string .= add_to_job_id_dependency_string(
                {
                    job_id_href       => $job_id_href,
                    case_id_chain_key => $case_id_chain_key,
                    chain_key         => $sample_id_chain_key,
                }
            );
        }

      PARALLEL_PROCESS:
        foreach my $parallel_processes_index (
            0 .. $max_parallel_processes_count_href->{$sample_id} )
        {

            # Create key
            my $sample_id_parallel_chain_key =
                $sample_id
              . $UNDERSCORE
              . q{parallel}
              . $UNDERSCORE
              . $path
              . $parallel_processes_index;

            ## If parallel job exists
            next PARALLEL_PROCESS
              if (
                not $job_id_href->{$case_id_chain_key}{$sample_id_parallel_chain_key} );

            ## Add to job_id string
            $job_ids_string .= add_to_job_id_dependency_string(
                {
                    job_id_href       => $job_id_href,
                    case_id_chain_key => $case_id_chain_key,
                    chain_key         => $sample_id_parallel_chain_key,
                }
            );
        }
    }
    return $job_ids_string;
}

sub add_parallel_job_ids_to_job_id_dependency_string {

## Function : Create job id string from the case parallel job id chain and path associated with case chain for SLURM submission using dependencies
## Returns  : $job_ids_string
## Arguments: $case_id_chain_key => Case id chain hash key
##          : $job_id_href       => The info on job ids hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;

    if ( defined $job_id_href->{$case_id_chain_key} ) {

      CASE_CHAIN_KEY:
        foreach my $chain_key ( keys %{ $job_id_href->{$case_id_chain_key} } ) {

            ## Check if chain_key actually is a parallel
            if ( $chain_key =~ /parallel/sxm ) {

                ## Add to job_id string
                $job_ids_string .= add_to_job_id_dependency_string(
                    {
                        case_id_chain_key => $case_id_chain_key,
                        chain_key         => $chain_key,
                        job_id_href       => $job_id_href,
                    }
                );
            }
        }
    }
    return $job_ids_string;
}

sub add_parallel_chains_job_ids_to_job_id_dependency_string {

## Function : Create job id string from the job id chain and path associated with parallel chains for SLURM submission using dependencies
## Returns  : $job_ids_string
## Arguments: $case_id             => Case id
##          : $job_id_href         => The info on job ids hash {REF}
##          : $parallel_chains_ref => Info on parallel chains array {REF}
##          : $sample_ids_ref      => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $job_id_href;
    my $parallel_chains_ref;
    my $sample_ids_ref;

    my $tmpl = {
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;

  PARALLEL_CHAINS:
    foreach my $parallel_chain ( @{$parallel_chains_ref} ) {

        my $case_id_parallel_chain_key = $case_id . $UNDERSCORE . $parallel_chain;
      SAMPLE_IDS:
        foreach my $sample_id ( @{$sample_ids_ref} ) {

            my $sample_id_parallel_chain_key = $sample_id . $UNDERSCORE . $parallel_chain;

            ## Add to job_id string
            $job_ids_string .= add_to_job_id_dependency_string(
                {
                    case_id_chain_key => $case_id_parallel_chain_key,
                    chain_key         => $sample_id_parallel_chain_key,
                    job_id_href       => $job_id_href,
                }
            );
        }

        ###  Family
        ## Add to job_id string
        $job_ids_string .= add_to_job_id_dependency_string(
            {
                case_id_chain_key => $case_id_parallel_chain_key,
                chain_key         => $case_id_parallel_chain_key,
                job_id_href       => $job_id_href,
            }
        );
    }
    return $job_ids_string;
}

sub add_job_id_dependency_tree {

## Function : Saves job_id to the correct hash array depending on chain type.
## Returns  :
## Arguments: $chain_key       => Arbitrary chain hash key
##          : $job_id_href     => Info on jobIds hash {REF}
##          : $job_id_returned => Job_id that was returned from submission

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_key;
    my $job_id_href;
    my $job_id_returned;

    my $tmpl = {
        chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        job_id_returned => {
            defined     => 1,
            required    => 1,
            store       => \$job_id_returned,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Push job_id to arbitrary chain key

    # Alias job_ids array for arbitrary chain key in job_id_href to push to
    my $chain_key_job_ids_ref =
      \@{ $job_id_href->{$chain_key}{$chain_key} };

    # Add to sample_id job dependency tree
    push @{$chain_key_job_ids_ref}, $job_id_returned;

    return;
}

sub add_parallel_job_id_to_sample_id_dependency_tree {

## Function : Saves job_id to the correct hash array in dependency tree hash depending on chain type
## Returns  :
## Arguments: $case_id_chain_key                 => Case id chain hash key
##          : $job_id_href                       => Info on job_ids hash {REF}
##          : $max_parallel_processes_count_href => Maximum number of parallel processes
##          : $path                              => Trunk or branch
##          : $sample_id                         => Sample ID
##          : $sample_id_chain_key                => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;
    my $max_parallel_processes_count_href;
    my $path;
    my $sample_id;
    my $sample_id_chain_key;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        max_parallel_processes_count_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$max_parallel_processes_count_href,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_id           => { store => \$sample_id, strict_type => 1, },
        sample_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id_chain_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $chain_key_type = q{parallel};
    my $parallel_jobs_chain_key;

    ## Push parallel job_ids
  PARALLEL_PROCESS:
    foreach my $parallel_processes_index (
        0 .. $max_parallel_processes_count_href->{$sample_id} )
    {

        # Set key
        $parallel_jobs_chain_key =
            $sample_id
          . $UNDERSCORE
          . $chain_key_type
          . $UNDERSCORE
          . $path
          . $parallel_processes_index;

        ## If parallel job_ids exists
        next PARALLEL_PROCESS
          if ( not exists $job_id_href->{$case_id_chain_key}{$parallel_jobs_chain_key} );

        # Alias job_ids array for sample_id in job_id_href to push to
        my $sample_id_job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$sample_id_chain_key} };

        # Alias parallel job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$parallel_jobs_chain_key} };

        push @{$sample_id_job_ids_ref}, @{$job_ids_ref};
    }
    return;
}

sub add_parallel_job_id_to_case_id_dependency_tree {

## Function : Saves job_id to the correct hash array in dependency tree hash depending on chain type.
## Returns  :
## Arguments: $case_id               => Case id
##          : $case_id_chain_key     => Case id chain hash key
##          : $job_id_href           => Info on jobIds hash {REF}
##          : $path                  => Trunk or branch
##          : $sbatch_script_tracker => Track the number of parallel processes (e.g. sbatch scripts for a module)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $case_id_chain_key;
    my $job_id_href;
    my $path;
    my $sbatch_script_tracker;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
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
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sbatch_script_tracker => {
            allow       => qr{ \A\d+\z }sxm,
            defined     => 1,
            required    => 1,
            store       => \$sbatch_script_tracker,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $chain_key_type = q{parallel};

    # Family parallel chainkey
    my $case_id_parallel_chain_key =
      $case_id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $sbatch_script_tracker;

    ## If parellel job_ids exists
    if ( exists $job_id_href->{$case_id_chain_key}{$case_id_parallel_chain_key} ) {

        # Alias job_ids array for sample_id in job_id_href to push to
        my $case_id_job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$case_id_chain_key} };

        # Alias parallel job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$case_id_parallel_chain_key} };

        push @{$case_id_job_ids_ref}, @{$job_ids_ref};
    }
    return;
}

sub add_parallel_job_id_to_parallel_dependency_tree {

## Function : Saves parallel job_id for case to the parallel case id hash array depending on chain type.
## Returns  :
## Arguments: $case_id_chain_key     => Case id chain hash key
##          : $id                    => Id (case or sample)
##          : $job_id_href           => Info on jobIds hash {REF}
##          : $job_id_returned       => Job_id that was returned from submission
##          : $path                  => Trunk or branch
##          : $sbatch_script_tracker => Track the number of parallel processes (e.g. sbatch scripts for a module)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $id;
    my $job_id_href;
    my $job_id_returned;
    my $path;
    my $sbatch_script_tracker;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        job_id_returned => {
            defined     => 1,
            required    => 1,
            store       => \$job_id_returned,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sbatch_script_tracker => {
            allow       => qr{ \A\d+\z }sxm,
            defined     => 1,
            required    => 1,
            store       => \$sbatch_script_tracker,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Push job to parallel job id

    my $chain_key_type = q{parallel};

    # Family parallel chainkey
    my $id_parallel_chain_key =
      $id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $sbatch_script_tracker;

    # Alias job_ids array for id in job_id_href to push to
    my $id_job_ids_ref =
      \@{ $job_id_href->{$case_id_chain_key}{$id_parallel_chain_key} };

    # Add to sample_id job dependency tree
    push @{$id_job_ids_ref}, $job_id_returned;

    return;
}

sub add_sample_id_parallel_job_id_to_case_id_dependency_tree {

## Function : Saves job_id to the correct hash array in dependency tree hash depending on chain type.
## Returns  :
## Arguments: $case_id_chain_key                 => Case id chain hash key
##          : $max_parallel_processes_count_href => Maximum number of parallel processes
##          : $job_id_href                       => Info on jobIds hash {REF}
##          : $path                              => Trunk or branch
##          : $sample_id                         => Sample ID

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;
    my $max_parallel_processes_count_href;
    my $path;
    my $sample_id;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        max_parallel_processes_count_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$max_parallel_processes_count_href,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_id => { store => \$sample_id, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not exists $max_parallel_processes_count_href->{$sample_id} );

    my $chain_key_type = q{parallel};
    my $parallel_jobs_chain_key;

    ## Push parallel job_ids
  PARALLEL_PROCESS:
    foreach my $parallel_processes_index (
        0 .. $max_parallel_processes_count_href->{$sample_id} )
    {

        # Set parallel sample key
        $parallel_jobs_chain_key =
            $sample_id
          . $UNDERSCORE
          . $chain_key_type
          . $UNDERSCORE
          . $path
          . $parallel_processes_index;

        ## If parallel job_ids exists
        next PARALLEL_PROCESS
          if ( not exists $job_id_href->{$case_id_chain_key}{$parallel_jobs_chain_key} );

        # Alias job_ids array for case_id in job_id_href to push to
        my $case_id_job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$case_id_chain_key} };

        # Alias parallel job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$parallel_jobs_chain_key} };

        push @{$case_id_job_ids_ref}, @{$job_ids_ref};
    }
    return;
}

sub add_sample_ids_parallel_job_id_to_case_id_dependency_tree {

## Function : Saves job_id to the correct hash array in dependency tree hash depending on chain type.
## Returns  :
## Arguments: $case_id_chain_key                 => Case id chain hash key
##          : $job_id_href                       => Info on jobIds hash {REF}
##          : $max_parallel_processes_count_href => Maximum number of parallel processes
##          : $path                              => Trunk or branch
##          : $sample_ids_ref                    => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;
    my $max_parallel_processes_count_href;
    my $path;
    my $sample_ids_ref;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        max_parallel_processes_count_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$max_parallel_processes_count_href,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        add_sample_id_parallel_job_id_to_case_id_dependency_tree(
            {
                case_id_chain_key                 => $case_id_chain_key,
                job_id_href                       => $job_id_href,
                max_parallel_processes_count_href => $max_parallel_processes_count_href,
                path                              => $path,
                sample_id                         => $sample_id,
            }
        );
    }
    return;
}

sub add_pan_job_id_to_sample_id_dependency_tree {

## Function : Saves pan (i.e job_ids that affect all chains) job_id to the the sample_id chain.
## Returns  :
## Arguments: $case_id_chain_key   => Case id chain hash key
##          : $job_id_href         => Info on jobIds hash {REF}
##          : $sample_id_chain_key => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;
    my $sample_id_chain_key;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        sample_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id_chain_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $pan_chain_key;

    ## Push pan jobs

    # Set key
    $pan_chain_key = $case_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    # If pan job ids exists
    if ( exists $job_id_href->{$case_id_chain_key}{$pan_chain_key} ) {

        # Alias pan job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$pan_chain_key} };

        # Alias job_ids array for sample_id in job_id_href to push to
        my $sample_id_job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$sample_id_chain_key} };

        # Add job_ids to sample_id job_ids
        push @{$sample_id_job_ids_ref}, @{$job_ids_ref};
    }
    return;
}

sub add_sample_job_id_to_sample_id_dependency_tree {

## Function : Saves sample job_id to the sample_id hash array depending on chain type.
## Returns  :
## Arguments: $case_id_chain_key   => Case id chain hash key
##          : $job_id_href         => Info on jobIds hash {REF}
##          : $job_id_returned     => Job_id that was returned from submission
##          : $sample_id_chain_key => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;
    my $job_id_returned;
    my $sample_id_chain_key;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        job_id_returned => {
            defined     => 1,
            required    => 1,
            store       => \$job_id_returned,
            strict_type => 1,
        },
        sample_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id_chain_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Push sample job

    # Alias job_ids array for sample_id in job_id_href to push to
    my $sample_id_job_ids_ref =
      \@{ $job_id_href->{$case_id_chain_key}{$sample_id_chain_key} };

    # Add to sample_id job dependency tree
    push @{$sample_id_job_ids_ref}, $job_id_returned;

    return;
}

sub add_sample_job_id_to_case_id_dependency_tree {

## Function : Saves sample_id job_id to the the case_id chain.
## Returns  :
## Arguments: $case_id_chain_key   => Case id chain hash key
##          : $job_id_href         => Info on jobIds hash {REF}
##          : $sample_id_chain_key => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;
    my $sample_id_chain_key;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        sample_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id_chain_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Push sample_id jobs

    if ( exists $job_id_href->{$case_id_chain_key}{$sample_id_chain_key} ) {

        # Alias case_pan job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$sample_id_chain_key} };

        ## Use $case_id_chain_key instead of $sample_id_chain_key
        # Alias job_ids array for case_id in job_id_href to push to
        my $case_id_job_ids_ref =
          \@{ $job_id_href->{$case_id_chain_key}{$case_id_chain_key} };

        # Add job_ids to case_id job_ids
        push @{$case_id_job_ids_ref}, @{$job_ids_ref};
    }
    return;
}

sub create_job_id_string_for_sample_id {

## Function : Create job id string from the job id chain and path associated with sample for SLURM submission using dependencies
## Returns  : $job_ids_string
## Arguments: $case_id               => Case id
##          : $case_id_chain_key     => Case id chain hash key
##          : $job_id_href           => The info on job ids hash {REF}
##          : $path                  => Trunk or branch
##          : $sample_id             => Sample id
##          : $sample_id_chain_key   => Sample ID chain hash key
##          : $sbatch_script_tracker => Track the number of parallel processes (e.g. sbatch scripts for a module)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $case_id_chain_key;
    my $job_id_href;
    my $path;
    my $sample_id;
    my $sample_id_chain_key;
    my $sbatch_script_tracker;

    my $tmpl = {
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        path      => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id_chain_key,
            strict_type => 1,
        },
        sbatch_script_tracker => {
            allow       => qr{ \A\d+\z }sxm,
            defined     => 1,
            store       => \$sbatch_script_tracker,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;
    my $path_main = q{MAIN};

    ## If MAIN path
    if ( $path eq $path_main ) {

        ## Add to job_id string
        $job_ids_string = add_to_job_id_dependency_string(
            {
                case_id_chain_key => $case_id_chain_key,
                chain_key         => $sample_id_chain_key,
                job_id_href       => $job_id_href,
            }
        );
    }
    elsif ( $path ne $path_main ) {

        ## NOT MAIN. Branch
        ## Check for any previous job_ids within path current PATH. Branch.

        my $case_id_chain_key_main   = $case_id . $UNDERSCORE . $path_main;
        my $sample_id_chain_key_main = $sample_id . $UNDERSCORE . $path_main;

        ## For sample parallel MAIN jobs
        my $sample_id_parallel_chain_key_main;

        # Inheritance from MAIN parallel jobs
        if ( defined $sbatch_script_tracker ) {

            $sample_id_parallel_chain_key_main =
                $sample_id
              . $UNDERSCORE
              . q{parallel}
              . $UNDERSCORE
              . $path_main
              . $sbatch_script_tracker;
        }
        ## Second or later in branch chain
        if ( $job_id_href->{$case_id_chain_key}{$sample_id_chain_key} ) {

            ## Add to job_id string
            $job_ids_string = add_to_job_id_dependency_string(
                {
                    case_id_chain_key => $case_id_chain_key,
                    chain_key         => $sample_id_chain_key,
                    job_id_href       => $job_id_href,
                }
            );
        }
        elsif ( $job_id_href->{$case_id_chain_key_main}{$sample_id_chain_key_main} ) {
            ## No previous job_ids with current path.
            ## Inherit from potential MAIN. Trunk

            ## Add to job_id string
            $job_ids_string = add_to_job_id_dependency_string(
                {
                    case_id_chain_key => $case_id_chain_key_main,
                    chain_key         => $sample_id_chain_key_main,
                    job_id_href       => $job_id_href,
                }
            );
        }
        elsif ( $sample_id_parallel_chain_key_main
            and
            $job_id_href->{$case_id_chain_key_main}{$sample_id_parallel_chain_key_main} )
        {
            ## No previous job_ids within MAIN path.
            ## Inherit from potential parallel jobs MAIN. Trunk

            ## Add to job_id string
            $job_ids_string = add_to_job_id_dependency_string(
                {
                    case_id_chain_key => $case_id_chain_key_main,
                    chain_key         => $sample_id_parallel_chain_key_main,
                    job_id_href       => $job_id_href,
                }
            );
        }
    }
    return $job_ids_string;
}

sub create_job_id_string_for_case_id {

## Function : Create job id string from the job id chain and path associated with case for SLURM submission using dependencies
## Returns  : $job_ids_string
##          : $case_id                           => Case id
##          : $case_id_chain_key                 => Case id chain hash key
##          : $job_id_href                       => The info on job ids hash {REF}
##          : $max_parallel_processes_count_href => Maximum number of parallel processes
##          : $parallel_chains_ref               => Info on parallel chains array {REF}
##          : $path                              => Trunk or branch
##          : $sample_ids_ref                    => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $case_id_chain_key;
    my $job_id_href;
    my $max_parallel_processes_count_href;
    my $parallel_chains_ref;
    my $path;
    my $sample_ids_ref;

    my $tmpl = {
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        max_parallel_processes_count_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$max_parallel_processes_count_href,
            strict_type => 1,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;
    my $path_main = q{MAIN};

    $job_ids_string = create_job_id_string_for_case_id_and_path(
        {
            case_id                           => $case_id,
            case_id_chain_key                 => $case_id_chain_key,
            job_id_href                       => $job_id_href,
            max_parallel_processes_count_href => $max_parallel_processes_count_href,
            path                              => $path,
            sample_ids_ref                    => $sample_ids_ref,
        }
    );

    if ( $path ne $path_main && !$job_ids_string ) {

        ## No previous sample_id job_ids for other path

        my $case_id_chain_key_main = $case_id . $UNDERSCORE . $path_main;

        $job_ids_string = create_job_id_string_for_case_id_and_path(
            {
                case_id                           => $case_id,
                case_id_chain_key                 => $case_id_chain_key_main,
                job_id_href                       => $job_id_href,
                max_parallel_processes_count_href => $max_parallel_processes_count_href,
                path                              => $path_main,
                sample_ids_ref                    => $sample_ids_ref,
            }
        );
    }
    if ( @{$parallel_chains_ref} ) {

        $job_ids_string .= add_parallel_chains_job_ids_to_job_id_dependency_string(
            {
                case_id             => $case_id,
                job_id_href         => $job_id_href,
                parallel_chains_ref => $parallel_chains_ref,
                sample_ids_ref      => $sample_ids_ref,
            }
        );

    }
    return $job_ids_string;
}

sub create_job_id_string_for_case_id_and_path {

## Function : Create job id string from the job id chain and main path for SLURM submission using dependencies
## Returns  : $job_ids_string
## Arguments: $case_id                           => Case id
##          : $case_id_chain_key                 => Case id chain hash key
##          : $job_id_href                       => The info on job ids hash {REF}
##          : $max_parallel_processes_count_href => Maximum number of parallel processes
##          : $path                              => Trunk or branch
##          : $sample_ids_ref                    => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $case_id_chain_key;
    my $job_id_href;
    my $max_parallel_processes_count_href;
    my $path;
    my $sample_ids_ref;

    my $tmpl = {
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        max_parallel_processes_count_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$max_parallel_processes_count_href,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;

    ## If previous case chain jobs exists, sample_ids have already been inherited
    if ( $job_id_href->{$case_id_chain_key}{$case_id_chain_key} ) {

        ## Add to job_id string
        $job_ids_string = add_to_job_id_dependency_string(
            {
                case_id_chain_key => $case_id_chain_key,
                chain_key         => $case_id_chain_key,
                job_id_href       => $job_id_href,
            }
        );
        return $job_ids_string;
    }

    ## First case_id in MAIN chain

    ## Add both parallel and sample_id job_id(s) from sample_id(s) chain_key
    $job_ids_string = add_sample_ids_job_ids_to_job_id_dependency_string(
        {
            case_id                           => $case_id,
            case_id_chain_key                 => $case_id_chain_key,
            job_id_href                       => $job_id_href,
            max_parallel_processes_count_href => $max_parallel_processes_count_href,
            path                              => $path,
            sample_ids_ref                    => $sample_ids_ref,
        }
    );
    return $job_ids_string;
}

sub clear_sample_id_parallel_job_id_dependency_tree {

## Function : Clear parallel sample job_ids in the sample_id chain
## Returns  :
## Arguments: $case_id_chain_key                 => Case id chain hash key
##          : $job_id_href                       => Info on jobIds hash {REF}
##          : $max_parallel_processes_count_href => Maximum number of parallel processes
##          : $path                              => Trunk or branch
##          : $sample_id                         => Sample ID

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;
    my $max_parallel_processes_count_href;
    my $path;
    my $sample_id;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        max_parallel_processes_count_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$max_parallel_processes_count_href,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
        sample_id => { store => \$sample_id, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Clear all latest parallel jobs within chain_key
  PARALLEL_PROCESS:
    foreach my $parallel_processes_index (
        0 .. $max_parallel_processes_count_href->{$sample_id} )
    {

        # Create key
        my $sample_id_parallel_chain_key =
            $sample_id
          . $UNDERSCORE
          . q{parallel}
          . $UNDERSCORE
          . $path
          . $parallel_processes_index;

        ## If parallel job
        if ( $job_id_href->{$case_id_chain_key}{$sample_id_parallel_chain_key} ) {

            # Clear latest parallel sample_id chain submission
            @{ $job_id_href->{$case_id_chain_key}{$sample_id_parallel_chain_key} } = ();
        }
    }
    return;
}

sub clear_pan_job_id_dependency_tree {

## Function : Clear pan (i.e job_ids that affect all sample chains) job_id to the pan chain.
## Returns  :
## Arguments: $case_id_chain_key   => Case id chain hash key
##          : $job_id_href         => Info on jobIds hash {REF}
##          : $sample_id_chain_key => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;
    my $sample_id_chain_key;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        sample_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id_chain_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $pan_chain_key;

    # Set key
    $pan_chain_key = $case_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    ## Clear latest case_id_sample_id chainkey
    @{ $job_id_href->{$case_id_chain_key}{$pan_chain_key} } = ();

    return;
}

sub clear_sample_id_job_id_dependency_tree {

## Function : Clear sample job_ids in the sample_id chain.
## Returns  :
## Arguments: $case_id_chain_key   => Case id chain hash key
##          : $job_id_href         => Info on jobIds hash {REF}
##          : $sample_id_chain_key => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;
    my $sample_id_chain_key;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        sample_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id_chain_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Clear latest sample_id chainkey
    @{ $job_id_href->{$case_id_chain_key}{$sample_id_chain_key} } = ();

    return;
}

sub clear_case_id_job_id_dependency_tree {

## Function : Clear case job_ids in the the case_id chain.
## Returns  :
## Arguments: $case_id_chain_key => Case id chain hash key
##          : $job_id_href       => Info on jobIds hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Clear latest case_id chainkey
    @{ $job_id_href->{$case_id_chain_key}{$case_id_chain_key} } = ();

    return;
}

sub clear_all_job_ids_within_chain_key_dependency_tree {

## Function : Clear all job_ids in the chain key.
## Returns  :
## Arguments: $case_id_chain_key => Case id chain hash key
##          : $job_id_href       => Info on jobIds hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $job_id_href;

    my $tmpl = {
        case_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Clear all jobs within chainkey
  CHAIN_KEYS:
    foreach my $chain_key ( keys %{ $job_id_href->{$case_id_chain_key} } ) {

        ## Clear all case_id/sample_id chain submission for path
        @{ $job_id_href->{$case_id_chain_key}{$chain_key} } =
          ();
    }
    return;
}

sub get_all_job_ids {

## Function : Get all job_ids
## Returns  : @job_ids
## Arguments: $job_id_href => Job id hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;

    my $tmpl = {
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not exists $job_id_href->{ALL} );

    return @{ $job_id_href->{ALL}{ALL} };
}

sub limit_job_id_string {

## Function : Limit number of job_ids in job_id chain
## Returns  :
## Arguments: $case_id_chain_key => Case id chain hash key
##          : $chain_key         => The current chain hash key
##          : $job_id_href       => The info on job ids hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id_chain_key;
    my $chain_key;
    my $job_id_href;

    my $tmpl = {
        case_id_chain_key => {
            default     => qw{ALL},
            store       => \$case_id_chain_key,
            strict_type => 1,
        },
        chain_key => {
            default     => qw{ALL},
            store       => \$chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Set maximum job_ids to track limit
    Readonly my $MAX_JOB_IDS_TO_TRACK => 100;

    # Alias job_id chain array
    my $job_ids_ref = $job_id_href->{$case_id_chain_key}{$chain_key};

    ## Keeps the job_id string dependency within reasonable limits
    if (   ( defined $job_ids_ref )
        && ( scalar @{$job_ids_ref} >= $MAX_JOB_IDS_TO_TRACK ) )
    {

        # Remove oldest job_ids
        @{$job_ids_ref} = @{$job_ids_ref}[ 0 .. $MAX_JOB_IDS_TO_TRACK ];
    }
    return;
}

sub print_wait {

## Function : Calculates when to print "wait" statement and prints "wait" to supplied filehandle when adequate.
## Returns  : $process_batches_count
## Arguments: $filehandle            => filehandle to print "wait" statment to
##          : $max_process_number    => The maximum number of processes to be use before printing "wait" statement
##          : $process_batches_count => Scales the number of $max_process_number processs used after each print "wait" statement
##          : $process_counter       => The number of started processes

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $max_process_number;
    my $process_batches_count;
    my $process_counter;

    my $tmpl = {
        filehandle         => { defined => 1, required => 1, store => \$filehandle, },
        max_process_number => {
            defined     => 1,
            required    => 1,
            store       => \$max_process_number,
            strict_type => 1,
        },
        process_batches_count => {
            defined     => 1,
            required    => 1,
            store       => \$process_batches_count,
            strict_type => 1,
        },
        process_counter => {
            defined     => 1,
            required    => 1,
            store       => \$process_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Bash qw{gnu_wait};

    # Using only nr of processs eq the maximum number of process scaled by the batch count
    if ( $process_counter == $process_batches_count * $max_process_number ) {

        # Print wait statement to filehandle
        gnu_wait( { filehandle => $filehandle, } );
        say {$filehandle} $NEWLINE;

# Increase the maximum number of processs allowed to be used since "wait" was just printed
        $process_batches_count = $process_batches_count + 1;
    }
    return $process_batches_count;
}

sub submit_recipe {

## Function : Submit recipe depending on submission profile
## Returns  :
## Arguments: $base_command                      => Profile base command
##          : $case_id                           => Case id
##          : $dependency_method                 => Dependency method
##          : $job_id_chain                      => Chain id
##          : $job_id_href                       => The info on job ids hash {REF}
##          : $job_dependency_type               => Job dependency type
##          : $log                               => Log object
##          : $max_parallel_processes_count_href => Maximum number of parallel processes
##          : $parallel_chains_ref               => Info on parallel chains array {REF}
##          : $recipe_file_path                  => Recipe file path
##          : $recipe_files_tracker              => Track the number of parallel processes (e.g. recipe scripts for a module)
##          : $job_reservation_name              => Allocate resources from named reservation
##          : $sample_id                         => Sample id
##          : $sample_ids_ref                    => Sample ids {REF}
##          : $submission_profile                => Submission profile

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $dependency_method;
    my $job_id_chain;
    my $job_id_href;
    my $job_dependency_type;
    my $log;
    my $max_parallel_processes_count_href;
    my $parallel_chains_ref;
    my $recipe_file_path;
    my $recipe_files_tracker;
    my $job_reservation_name;
    my $sample_id;
    my $sample_ids_ref;

    ## Default(s)
    my $base_command;
    my $submission_profile;

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
            defined     => 1,
            required    => 1,
            store       => \$dependency_method,
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
        job_dependency_type => {
            store       => \$job_dependency_type,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        max_parallel_processes_count_href => {
            default     => {},
            store       => \$max_parallel_processes_count_href,
            strict_type => 1,
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
        recipe_files_tracker => {
            allow       => qr{ \A\d+\z }sxm,
            store       => \$recipe_files_tracker,
            strict_type => 1,
        },
        job_reservation_name => {
            store       => \$job_reservation_name,
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
        submission_profile => {
            allow       => [qw{ slurm }],
            default     => q{slurm},
            store       => \$submission_profile,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Slurm_processes qw{ submit_slurm_recipe };

    my %is_manager = ( slurm => \&submit_slurm_recipe, );

    $is_manager{$submission_profile}->(
        {
            base_command                      => $base_command,
            case_id                           => $case_id,
            dependency_method                 => $dependency_method,
            job_dependency_type               => $job_dependency_type,
            job_id_chain                      => $job_id_chain,
            job_id_href                       => $job_id_href,
            log                               => $log,
            max_parallel_processes_count_href => $max_parallel_processes_count_href,
            parallel_chains_ref               => $parallel_chains_ref,
            recipe_file_path                  => $recipe_file_path,
            recipe_files_tracker              => $recipe_files_tracker,
            reservation_name                  => $job_reservation_name,
            sample_id                         => $sample_id,
            sample_ids_ref                    => $sample_ids_ref,
        }
    );
    return 1;
}

sub write_job_ids_to_file {

## Function : Write all job_ids to file
## Returns  :
## Arguments: $case_id         => Case id
##          : $date_time_stamp => The date and time
##          : $log_file        => Log file
##          : $job_id_href     => Job id hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $date_time_stamp;
    my $log_file;
    my $job_id_href;

    my $tmpl = {
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        date_time_stamp => {
            defined     => 1,
            required    => 1,
            store       => \$date_time_stamp,
            strict_type => 1,
        },
        log_file => {
            defined     => 1,
            required    => 1,
            store       => \$log_file,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Write qw{ write_to_file };

    ## Write job_ids file
    return if ( not keys %{$job_id_href} );

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $log_dir      = dirname($log_file);
    my $job_ids_file = catfile( $log_dir,
        q{slurm_job_ids} . $UNDERSCORE . $date_time_stamp . $DOT . q{yaml} );

    ## Remove all undef elements from return array of all job_ids
    my @job_ids =
      grep { defined } get_all_job_ids( { job_id_href => $job_id_href, } );

    ## Writes a YAML hash to file
    my %out_job_id = ( $case_id => [@job_ids], );
    write_to_file(
        {
            data_href => \%out_job_id,
            format    => q{yaml},
            path      => $job_ids_file,
        }
    );
    $log->info( q{Wrote: } . $job_ids_file );
    return 1;
}

1;
