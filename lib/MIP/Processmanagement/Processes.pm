package MIP::Processmanagement::Processes;

use Carp;
use charnames qw{ :full :short };
use FindBin qw{$Bin};    # Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{check allow last_error};
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

BEGIN {
    use base qw (Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      add_to_job_id_dependency_string
      add_sample_ids_job_ids_to_job_id_dependency_string
      add_parallel_job_ids_to_job_id_dependency_string
      add_parallel_chains_job_ids_to_job_id_dependency_string
      add_job_id_dependency_tree
      add_parallel_job_id_to_sample_id_dependency_tree
      add_parallel_job_id_to_family_id_dependency_tree
      add_parallel_job_id_to_parallel_dependency_tree
      add_sample_id_parallel_job_id_to_family_id_dependency_tree
      add_sample_ids_parallel_job_id_to_family_id_dependency_tree
      add_pan_job_id_to_sample_id_dependency_tree
      add_sample_job_id_to_sample_id_dependency_tree
      add_sample_job_id_to_family_id_dependency_tree
      create_job_id_string_for_sample_id
      create_job_id_string_for_family_id
      clear_sample_id_parallel_job_id_dependency_tree
      clear_pan_job_id_dependency_tree
      clear_sample_id_job_id_dependency_tree
      clear_family_id_job_id_dependency_tree
      clear_all_job_ids_within_chain_key_dependency_tree
      limit_job_id_string
      print_wait};
}

##Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};
Readonly my $EMPTY_STR  => q{};
Readonly my $COLON      => q{:};

sub add_to_job_id_dependency_string {

##add_to_job_id_dependency_string

##Function : Adds all previous job_ids per family_chain_key and chain_key to job_ids dependency string, which is used to set the dependency in SLURM.
##Returns  : "$job_ids"
##Arguments: $job_id_href, $family_id_chain_key, $chain_key
##         : $job_id_href         => Info on jobIds hash {REF}
##         : $family_id_chain_key => Family ID chain hash key
##         : $chain_key           => The current chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;
    my $chain_key;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$chain_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Job_id string to submit to workload manager
    my $job_ids_string = $EMPTY_STR;

    # Alias hash to get job_id dependencies from
    my $chain_job_id_href = $job_id_href->{$family_id_chain_key}{$chain_key};

    if ($chain_job_id_href) {

      JOB_IDS:
        while ( my ( $job_index, $job_id ) = each @{$chain_job_id_href} ) {

            # Only for defined job_ids
            if ( defined $job_id ) {

                ## Only 1 previous job_id
                if ( ( !$job_index ) && ( scalar @{$chain_job_id_href} == 1 ) )
                {

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

##add_sample_ids_job_ids_to_job_id_dependency_string

##Function : Create job id string from sample_ids job id chain and path for SLURM submission using dependencies
##Returns  : "$job_ids_string"
##Arguments: $job_id_href, $infile_lane_prefix_href, $sample_ids_ref, $family_id, $family_id_chain_key, $path
##         : $job_id_href             => The info on job ids hash {REF}
##         : $infile_lane_prefix_href => The infile(s) without the ".ending" {REF}
##         : $sample_ids_ref          => Sample ids {REF}
##         : $family_id               => Family id
##         : $family_id_chain_key     => Family ID chain hash key
##         : $path                    => Trunk or branch

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $infile_lane_prefix_href;
    my $sample_ids_ref;
    my $family_id_chain_key;
    my $family_id;
    my $path;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
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
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;

    ## Add all previous jobId(s) from sample_id chainkey(s)
  SAMPLE_IDS:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

        if ( $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} ) {

            ## Add to job_id string
            $job_ids_string .= add_to_job_id_dependency_string(
                {
                    job_id_href         => $job_id_href,
                    family_id_chain_key => $family_id_chain_key,
                    chain_key           => $sample_id_chain_key,
                }
            );
        }

      INFILES:
        while ( my ($infile_index) =
            each @{ $infile_lane_prefix_href->{$sample_id} } )
        {

            # Create key
            my $sample_id_parallel_chain_key =
                $sample_id
              . $UNDERSCORE
              . q{parallel}
              . $UNDERSCORE
              . $path
              . $infile_index;

            ## If parallel job exists
            if ( $job_id_href->{$family_id_chain_key}
                {$sample_id_parallel_chain_key} )
            {

                ## Add to job_id string
                $job_ids_string .= add_to_job_id_dependency_string(
                    {
                        job_id_href         => $job_id_href,
                        family_id_chain_key => $family_id_chain_key,
                        chain_key           => $sample_id_parallel_chain_key,
                    }
                );
            }
        }
    }
    return $job_ids_string;
}

sub add_parallel_job_ids_to_job_id_dependency_string {

##add_parallel_job_ids_to_job_id_dependency_string

##Function : Create job id string from the family parallel job id chain and path associated with family chain for SLURM submission using dependencies
##Returns  : "$job_ids_string"
##Arguments: $job_id_href, $family_id_chain_key
##         : $job_id_href         => The info on job ids hash {REF}
##         : $family_id_chain_key => Family ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;

    if ( defined $job_id_href->{$family_id_chain_key} ) {

        foreach my $chain_key ( keys %{ $job_id_href->{$family_id_chain_key} } )
        {

            ## Check if chain_key actually is a parallel
            if ( $chain_key =~ /parallel/ ) {

                ## Add to job_id string
                $job_ids_string .= add_to_job_id_dependency_string(
                    {
                        job_id_href         => $job_id_href,
                        family_id_chain_key => $family_id_chain_key,
                        chain_key           => $chain_key,
                    }
                );
            }
        }
    }
    return $job_ids_string;
}

sub add_parallel_chains_job_ids_to_job_id_dependency_string {

##add_parallel_chains_job_ids_to_job_id_dependency_string

##Function : Create job id string from the job id chain and path associated with parallel chains for SLURM submission using dependencies
##Returns  : "$job_ids_string"
##Arguments: $job_id_href, $sample_ids_ref, $parallel_chains_ref, $family_id
##         : $job_id_href         => The info on job ids hash {REF}
##         : $sample_ids_ref      => Sample ids {REF}
##         : $parallel_chains_ref => Info on parallel chains array {REF}
##         : $family_id           => Family id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $sample_ids_ref;
    my $parallel_chains_ref;
    my $family_id;

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
        parallel_chains_ref => {
            default     => [],
            strict_type => 1,
            store       => \$parallel_chains_ref
        },
        family_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;

  PARALLEL_CHAINS:
    foreach my $parallel_chain ( @{$parallel_chains_ref} ) {

        my $family_id_parallel_chain_key =
          $family_id . $UNDERSCORE . $parallel_chain;

      SAMPLE_IDS:
        foreach my $sample_id ( @{$sample_ids_ref} ) {

            my $sample_id_parallel_chain_key =
              $sample_id . $UNDERSCORE . $parallel_chain;

            ## Add to job_id string
            $job_ids_string .= add_to_job_id_dependency_string(
                {
                    job_id_href         => $job_id_href,
                    family_id_chain_key => $family_id_parallel_chain_key,
                    chain_key           => $sample_id_parallel_chain_key,
                }
            );
        }

        ###  Family
        ## Add to job_id string
        $job_ids_string .= add_to_job_id_dependency_string(
            {
                job_id_href         => $job_id_href,
                family_id_chain_key => $family_id_parallel_chain_key,
                chain_key           => $family_id_parallel_chain_key,
            }
        );
    }
    return $job_ids_string;
}

sub add_job_id_dependency_tree {

##add_job_id_dependency_tree

##Function : Saves job_id to the correct hash array depending on chain type.
##Returns  : ""
##Arguments: $job_id_href, $chain_key, $job_id_returned
##         : $job_id_href     => Info on jobIds hash {REF}
##         : $chain_key       => Arbitrary chain hash key
##         : $job_id_returned => Job_id that was returned from submission

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $chain_key;
    my $job_id_returned;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$chain_key
        },
        job_id_returned => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$job_id_returned
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

##add_parallel_job_id_to_sample_id_dependency_tree

##Function : Saves job_id to the correct hash array in dependency tree hash depending on chain type.
##Returns  : ""
##Arguments: $job_id_href, $infile_lane_prefix_href, $family_id_chain_key, $sample_id_chain_key, $sample_id, $path
##         : $job_id_href             => Info on jobIds hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $family_id_chain_key     => Family ID chain hash key
##         : $sample_id_chain_key     => Sample ID chain hash key
##         : $sample_id               => Sample ID
##         : $path                    => Trunk or branch

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $family_id_chain_key;
    my $sample_id_chain_key;
    my $sample_id;
    my $path;

    my $tmpl = {
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        sample_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id_chain_key
        },
        sample_id => { strict_type => 1, store => \$sample_id },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $chain_key_type = q{parallel};
    my $parallel_jobs_chain_key;

    ## Push parallel job_ids
  INFILES:
    while ( my ($infile_index) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        # Set key
        $parallel_jobs_chain_key =
            $sample_id
          . $UNDERSCORE
          . $chain_key_type
          . $UNDERSCORE
          . $path
          . $infile_index;

        ## If parellel job_ids exists
        if (
            exists $job_id_href->{$family_id_chain_key}
            {$parallel_jobs_chain_key} )
        {

            # Alias job_ids array for sample_id in job_id_href to push to
            my $sample_id_job_ids_ref =
              \@{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} };

            # Alias parallel job_ids array to push from
            my $job_ids_ref =
              \@{ $job_id_href->{$family_id_chain_key}{$parallel_jobs_chain_key}
              };

            push @{$sample_id_job_ids_ref}, @{$job_ids_ref};
        }
    }
    return;
}

sub add_parallel_job_id_to_family_id_dependency_tree {

##add_parallel_job_id_to_family_id_dependency_tree

##Function : Saves job_id to the correct hash array in dependency tree hash depending on chain type.
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key, $family_id, $path, $sbatch_script_tracker
##         : $job_id_href           => Info on jobIds hash {REF}
##         : $family_id_chain_key   => Family ID chain hash key
##         : $family_id             => Family id
##         : $path                  => Trunk or branch
##         : $sbatch_script_tracker => Track the number of parallel processes (e.g. sbatch scripts for a module)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;
    my $family_id;
    my $path;
    my $sbatch_script_tracker;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        family_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
        sbatch_script_tracker => {
            required    => 1,
            defined     => 1,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$sbatch_script_tracker
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $chain_key_type = q{parallel};

    # Family parallel chainkey
    my $family_id_parallel_chain_key =
        $family_id
      . $UNDERSCORE
      . q{parallel}
      . $UNDERSCORE
      . $path
      . $sbatch_script_tracker;

    ## If parellel job_ids exists
    if (
        exists $job_id_href->{$family_id_chain_key}
        {$family_id_parallel_chain_key} )
    {

        # Alias job_ids array for sample_id in job_id_href to push to
        my $family_id_job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}{$family_id_chain_key} };

        # Alias parallel job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}
              {$family_id_parallel_chain_key} };

        push @{$family_id_job_ids_ref}, @{$job_ids_ref};
    }
    return;
}

sub add_parallel_job_id_to_parallel_dependency_tree {

##add_parallel_job_id_to_parallel_dependency_tree

##Function : Saves parallel job_id for family to the parallel family id hash array depending on chain type.
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key, $id, $path, $sbatch_script_tracker, $job_id_returned
##         : $job_id_href           => Info on jobIds hash {REF}
##         : $family_id_chain_key   => Family ID chain hash key
##         : $id                    => Id (family or sample)
##         : $path                  => Trunk or branch
##         : $sbatch_script_tracker => Track the number of parallel processes (e.g. sbatch scripts for a module)
##         : $job_id_returned       => Job_id that was returned from submission

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;
    my $id;
    my $path;
    my $sbatch_script_tracker;
    my $job_id_returned;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$id
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
        sbatch_script_tracker => {
            required    => 1,
            defined     => 1,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$sbatch_script_tracker
        },
        job_id_returned => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$job_id_returned
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Push job to parallel job id

    my $chain_key_type = q{parallel};

    # Family parallel chainkey
    my $id_parallel_chain_key =
        $id
      . $UNDERSCORE
      . q{parallel}
      . $UNDERSCORE
      . $path
      . $sbatch_script_tracker;

    # Alias job_ids array for id in job_id_href to push to
    my $id_job_ids_ref =
      \@{ $job_id_href->{$family_id_chain_key}{$id_parallel_chain_key} };

    # Add to sample_id job dependency tree
    push @{$id_job_ids_ref}, $job_id_returned;

    return;
}

sub add_sample_id_parallel_job_id_to_family_id_dependency_tree {

##add_sample_id_parallel_job_id_to_sample_id_dependency_tree

##Function : Saves job_id to the correct hash array in dependency tree hash depending on chain type.
##Returns  : ""
##Arguments: $job_id_href, $infile_lane_prefix_href, $family_id_chain_key, $sample_id, $path
##         : $job_id_href             => Info on jobIds hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $family_id_chain_key     => Family ID chain hash key
##         : $sample_id               => Sample ID
##         : $path                    => Trunk or branch

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $family_id_chain_key;
    my $sample_id_chain_key;
    my $sample_id;
    my $path;

    my $tmpl = {
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        sample_id => { strict_type => 1, store => \$sample_id },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $chain_key_type = q{parallel};
    my $parallel_jobs_chain_key;

    if ( exists $infile_lane_prefix_href->{$sample_id} ) {

        ## Push parallel job_ids
      INFILES:
        while ( my ($infile_index) =
            each @{ $infile_lane_prefix_href->{$sample_id} } )
        {

            # Set parallel sample key
            $parallel_jobs_chain_key =
                $sample_id
              . $UNDERSCORE
              . $chain_key_type
              . $UNDERSCORE
              . $path
              . $infile_index;

            ## If parellel job_ids exists
            if (
                exists $job_id_href->{$family_id_chain_key}
                {$parallel_jobs_chain_key} )
            {

                # Alias job_ids array for family_id in job_id_href to push to
                my $family_id_job_ids_ref =
                  \@{ $job_id_href->{$family_id_chain_key}{$family_id_chain_key}
                  };

                # Alias parallel job_ids array to push from
                my $job_ids_ref =
                  \@{ $job_id_href->{$family_id_chain_key}
                      {$parallel_jobs_chain_key} };

                push @{$family_id_job_ids_ref}, @{$job_ids_ref};
            }
        }
    }
    return;
}

sub add_sample_ids_parallel_job_id_to_family_id_dependency_tree {

##add_sample_ids_parallel_job_id_to_family_id_dependency_tree

##Function : Saves job_id to the correct hash array in dependency tree hash depending on chain type.
##Returns  : ""
##Arguments: $job_id_href, $infile_lane_prefix_href, $sample_ids_ref, $family_id_chain_key, $sample_id_chain_key, $sample_id, $path
##         : $job_id_href             => Info on jobIds hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $sample_ids_ref          => Sample ids {REF}
##         : $family_id_chain_key     => Family ID chain hash key
##         : $sample_id_chain_key     => Sample ID chain hash key
##         : $sample_id               => Sample ID
##         : $path                    => Trunk or branch

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_ids_ref;
    my $family_id_chain_key;
    my $sample_id_chain_key;
    my $sample_id;
    my $path;

    my $tmpl = {
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
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
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  SAMPLE_IDS:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        add_sample_id_parallel_job_id_to_family_id_dependency_tree(
            {
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                family_id_chain_key     => $family_id_chain_key,
                sample_id               => $sample_id,
                path                    => $path,
            }
        );
    }
    return;
}

sub add_pan_job_id_to_sample_id_dependency_tree {

##add_pan_job_id_to_sample_id_dependency_tree

##Function : Saves pan (i.e job_ids that affect all chains) job_id to the the sample_id chain.
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key, $sample_id_chain_key
##         : $job_id_href         => Info on jobIds hash {REF}
##         : $family_id_chain_key => Family ID chain hash key
##         : $sample_id_chain_key => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;
    my $sample_id_chain_key;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        sample_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id_chain_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $pan_chain_key;

    ## Push pan jobs

    # Set key
    $pan_chain_key = $family_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    # If pan job ids exists
    if ( exists $job_id_href->{$family_id_chain_key}{$pan_chain_key} ) {

        # Alias pan job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}{$pan_chain_key} };

        # Alias job_ids array for sample_id in job_id_href to push to
        my $sample_id_job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} };

        # Add job_ids to sample_id job_ids
        push @{$sample_id_job_ids_ref}, @{$job_ids_ref};
    }
    return;
}

sub add_sample_job_id_to_sample_id_dependency_tree {

##add_sample_job_id_to_sample_id_dependency_tree

##Function : Saves sample job_id to the sample_id hash array depending on chain type.
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key, $sample_id_chain_key, $job_id_returned
##         : $job_id_href         => Info on jobIds hash {REF}
##         : $family_id_chain_key => Family ID chain hash key
##         : $sample_id_chain_key => Sample ID chain hash key
##         : $job_id_returned     => Job_id that was returned from submission

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;
    my $sample_id_chain_key;
    my $job_id_returned;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        sample_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id_chain_key
        },
        job_id_returned => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$job_id_returned
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Push sample job

    # Alias job_ids array for sample_id in job_id_href to push to
    my $sample_id_job_ids_ref =
      \@{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} };

    # Add to sample_id job dependency tree
    push @{$sample_id_job_ids_ref}, $job_id_returned;

    return;
}

sub add_sample_job_id_to_family_id_dependency_tree {

##add_sample_job_id_to_family_id_dependency_tree

##Function : Saves sample_id job_id to the the family_id chain.
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key, $sample_id_chain_key
##         : $job_id_href         => Info on jobIds hash {REF}
##         : $family_id_chain_key => Family ID chain hash key
##         : $sample_id_chain_key => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;
    my $sample_id_chain_key;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        sample_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id_chain_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Push sample_id jobs

    if ( exists $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} ) {

        # Alias family_pan job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} };

        ## Use $family_id_chain_key instead of $sample_id_chain_key
        # Alias job_ids array for family_id in job_id_href to push to
        my $family_id_job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}{$family_id_chain_key} };

        # Add job_ids to family_id job_ids
        push @{$family_id_job_ids_ref}, @{$job_ids_ref};
    }
    return;
}

sub create_job_id_string_for_sample_id {

## Function : Create job id string from the job id chain and path associated with sample for SLURM submission using dependencies
## Returns  : $job_ids_string
## Arguments: $family_id             => Family id
##          : $family_id_chain_key   => Family ID chain hash key
##          : $job_id_href           => The info on job ids hash {REF}
##          : $path                  => Trunk or branch
##          : $sample_id             => Sample id
##          : $sample_id_chain_key   => Sample ID chain hash key
##          : $sbatch_script_tracker => Track the number of parallel processes (e.g. sbatch scripts for a module)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id;
    my $family_id_chain_key;
    my $path;
    my $sample_id;
    my $sample_id_chain_key;
    my $sbatch_script_tracker;

    my $tmpl = {
        family_id => {
            defined     => 1,
            required    => 1,
            store       => \$family_id,
            strict_type => 1,
        },
        family_id_chain_key => {
            defined     => 1,
            required    => 1,
            store       => \$family_id_chain_key,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        path =>
          { defined => 1, required => 1, store => \$path, strict_type => 1, },
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
            allow       => qr/^\d+$/,
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
                job_id_href         => $job_id_href,
                family_id_chain_key => $family_id_chain_key,
                chain_key           => $sample_id_chain_key,
            }
        );
    }
    elsif ( $path ne $path_main ) {
        ## NOT MAIN. Branch
        ## Check for any previous job_ids within path current PATH. Branch.

        my $sample_id_chain_key_main = $sample_id . $UNDERSCORE . $path_main;
        my $family_id_chain_key_main = $family_id . $UNDERSCORE . $path_main;
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
        elsif (
            $job_id_href->{$family_id_chain_key_main}{$sample_id_chain_key_main}
          )
        {
            ## No previous job_ids with current path.
            ## Inherit from potential MAIN. Trunk

            ## Add to job_id string
            $job_ids_string = add_to_job_id_dependency_string(
                {
                    job_id_href         => $job_id_href,
                    family_id_chain_key => $family_id_chain_key_main,
                    chain_key           => $sample_id_chain_key_main,
                }
            );
        }
        elsif ( $job_id_href->{$family_id_chain_key_main}
            {$sample_id_parallel_chain_key_main} )
        {
            ## No previous job_ids within MAIN path.
            ## Inherit from potential parallel jobs MAIN. Trunk

            ## Add to job_id string
            $job_ids_string = add_to_job_id_dependency_string(
                {
                    job_id_href         => $job_id_href,
                    family_id_chain_key => $family_id_chain_key_main,
                    chain_key           => $sample_id_parallel_chain_key_main,
                }
            );
        }
    }
    return $job_ids_string;
}

sub create_job_id_string_for_family_id {

##create_job_id_string_for_family_id

##Function : Create job id string from the job id chain and path associated with family for SLURM submission using dependencies
##Returns  : "$job_ids_string"
##Arguments: $job_id_href, $infile_lane_prefix_href, $sample_ids_ref, $parallel_chains_ref, $family_id, $family_id_chain_key, $path
##         : $job_id_href             => The info on job ids hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $sample_ids_ref          => Sample ids {REF}
##         : $parallel_chains_ref     => Info on parallel chains array {REF}
##         : $family_id               => Family id
##         : $family_id_chain_key     => Family ID chain hash key
##         : $path                    => Trunk or branch

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $infile_lane_prefix_href;
    my $sample_ids_ref;
    my $parallel_chains_ref;
    my $family_id_chain_key;
    my $family_id;
    my $path;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref
        },
        parallel_chains_ref => {
            default     => [],
            strict_type => 1,
            store       => \$parallel_chains_ref
        },
        family_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;
    my $path_main = q{MAIN};

    $job_ids_string = create_job_id_string_for_family_id_and_path(
        {
            job_id_href             => $job_id_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            sample_ids_ref          => $sample_ids_ref,
            family_id               => $family_id,
            family_id_chain_key     => $family_id_chain_key,
            path                    => $path,
        }
    );

    if ( $path ne $path_main && !$job_ids_string ) {
        ## No previous sample_id job_ids for other path

        my $family_id_chain_key_main = $family_id . $UNDERSCORE . $path_main;

        $job_ids_string = create_job_id_string_for_family_id_and_path(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref          => $sample_ids_ref,
                family_id               => $family_id,
                family_id_chain_key     => $family_id_chain_key_main,
                path                    => $path_main,
            }
        );
    }
    if ( @{$parallel_chains_ref} ) {

        $job_ids_string .=
          add_parallel_chains_job_ids_to_job_id_dependency_string(
            {
                job_id_href         => $job_id_href,
                sample_ids_ref      => $sample_ids_ref,
                parallel_chains_ref => $parallel_chains_ref,
                family_id           => $family_id,
            }
          );

    }
    return $job_ids_string;
}

sub create_job_id_string_for_family_id_and_path {

##create_job_id_string_for_family_id_and_path

##Function : Create job id string from the job id chain and main path for SLURM submission using dependencies
##Returns  : "$job_ids_string"
##Arguments: $job_id_href, $infile_lane_prefix_href, $sample_ids_ref, $family_id, $family_id_chain_key, $path
##         : $job_id_href             => The info on job ids hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $sample_ids_ref          => Sample ids {REF}
##         : $family_id               => Family id
##         : $family_id_chain_key     => Family ID chain hash key
##         : $path                    => Trunk or branch

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $infile_lane_prefix_href;
    my $sample_ids_ref;
    my $family_id_chain_key;
    my $family_id;
    my $path;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
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
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $job_ids_string;

    ## If previous family chain jobs exists, sample_ids have already been inherited
    if ( $job_id_href->{$family_id_chain_key}{$family_id_chain_key} ) {

        ## Add to job_id string
        $job_ids_string = add_to_job_id_dependency_string(
            {
                job_id_href         => $job_id_href,
                family_id_chain_key => $family_id_chain_key,
                chain_key           => $family_id_chain_key,
            }
        );
    }
    else {
        ## First family_id in MAIN chain

        ## Add both parallel and sample_id jobId(s) from sample_id(s) chainkey
        $job_ids_string = add_sample_ids_job_ids_to_job_id_dependency_string(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref          => $sample_ids_ref,
                family_id_chain_key     => $family_id_chain_key,
                family_id               => $family_id,
                path                    => $path,
            }
        );
    }
    return $job_ids_string;
}

sub clear_sample_id_parallel_job_id_dependency_tree {

##clear_sample_id_parallel_job_id_dependency_tree

##Function : Clear parallel sample job_ids in the sample_id chain.
##Returns  : ""
##Arguments: $job_id_href, $infile_lane_prefix_href, $family_id_chain_key, $sample_id, $path
##         : $job_id_href             => Info on jobIds hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $family_id_chain_key     => Family ID chain hash key
##         : $sample_id               => Sample ID
##         : $path                    => Trunk or branch

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $infile_lane_prefix_href;
    my $family_id_chain_key;
    my $sample_id;
    my $path;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        sample_id => { strict_type => 1, store => \$sample_id },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Clear all latest parallel jobs within chainkey
  INFILES:
    while ( my ($infile_index) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        # Create key
        my $sample_id_parallel_chain_key =
            $sample_id
          . $UNDERSCORE
          . q{parallel}
          . $UNDERSCORE
          . $path
          . $infile_index;

        ## If parallel job exists
        if ( $job_id_href->{$family_id_chain_key}{$sample_id_parallel_chain_key}
          )
        {

            # Clear latest parallel sample_id chain submission
            @{ $job_id_href->{$family_id_chain_key}
                  {$sample_id_parallel_chain_key} } = ();
        }
    }
    return;
}

sub clear_pan_job_id_dependency_tree {

##clear_pan_job_id_dependency_tree

##Function : Clear pan (i.e job_ids that affect all sample chains) job_id to the pan chain.
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key, $sample_id_chain_key
##         : $job_id_href         => Info on jobIds hash {REF}
##         : $family_id_chain_key => Family ID chain hash key
##         : $sample_id_chain_key => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;
    my $sample_id_chain_key;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        sample_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id_chain_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $pan_chain_key;

    # Set key
    $pan_chain_key = $family_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    ## Clear latest family_id_sample_id chainkey
    @{ $job_id_href->{$family_id_chain_key}{$pan_chain_key} } = ();

    return;
}

sub clear_sample_id_job_id_dependency_tree {

##clear_sample_id_job_id_dependency_tree

##Function : Clear sample job_ids in the sample_id chain.
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key, $sample_id_chain_key
##         : $job_id_href         => Info on jobIds hash {REF}
##         : $family_id_chain_key => Family ID chain hash key
##         : $sample_id_chain_key => Sample ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;
    my $sample_id_chain_key;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        sample_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id_chain_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Clear latest sample_id chainkey
    @{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} } = ();

    return;
}

sub clear_family_id_job_id_dependency_tree {

##clear_family_id_job_id_dependency_tree

##Function : Clear family job_ids in the the family_id chain.
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key
##         : $job_id_href         => Info on jobIds hash {REF}
##         : $family_id_chain_key => Family ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Clear latest family_id chainkey
    @{ $job_id_href->{$family_id_chain_key}{$family_id_chain_key} } = ();

    return;
}

sub clear_all_job_ids_within_chain_key_dependency_tree {

##clear_all_job_ids_within_chain_key_dependency_tree

##Function : Clear all job_ids in the chain key.
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key
##         : $job_id_href         => Info on jobIds hash {REF}
##         : $family_id_chain_key => Family ID chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_id_chain_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Clear all jobs within chainkey
  CHAIN_KEYS:
    foreach my $chain_key ( keys %{ $job_id_href->{$family_id_chain_key} } ) {

        ## Clear all family_id/sample_id chain submission for path
        @{ $job_id_href->{$family_id_chain_key}{$chain_key} } =
          ();
    }
    return;
}

sub limit_job_id_string {

##limit_job_id_string

##Function : Limit number of job_ids in job_id chain
##Returns  : ""
##Arguments: $job_id_href, $family_id_chain_key, $chain_key
##         : $job_id_href         => The info on job ids hash {REF}
##         : $family_id_chain_key => Family ID chain hash key
##         : $chain_key           => The current chain hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $family_id_chain_key;
    my $chain_key;

    my $tmpl = {
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        family_id_chain_key => {
            default     => qw{ALL},
            strict_type => 1,
            store       => \$family_id_chain_key
        },
        chain_key => {
            default     => qw{ALL},
            strict_type => 1,
            store       => \$chain_key
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Set maximum job_ids to track limit
    Readonly my $MAX_JOB_IDS_TO_TRACK => 100;

    # Alias job_id chain array
    my $job_ids_ref = $job_id_href->{$family_id_chain_key}{$chain_key};

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

##print_wait

##Function : Calculates when to print "wait" statement and prints "wait" to supplied FILEHANDLE when adequate.
##Returns  : "$process_batches_count"
##Arguments: $process_counter, $max_process_number, $process_batches_count, $FILEHANDLE
##         : $process_counter       => The number of started processes
##         : $max_process_number    => The maximum number of processes to be use before printing "wait" statement
##         : $process_batches_count => Scales the number of $max_process_number processs used after each print "wait" statement
##         : $FILEHANDLE            => FILEHANDLE to print "wait" statment to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $process_counter;
    my $max_process_number;
    my $process_batches_count;
    my $FILEHANDLE;

    my $tmpl = {
        process_counter => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$process_counter
        },
        max_process_number => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$max_process_number
        },
        process_batches_count => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$process_batches_count
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Bash qw{gnu_wait};

# Using only nr of processs eq the maximum number of process scaled by the batch count
    if ( $process_counter == $process_batches_count * $max_process_number ) {

        # Print wait statement to filehandle
        gnu_wait( { FILEHANDLE => $FILEHANDLE, } );
        say {$FILEHANDLE} $NEWLINE;

# Increase the maximum number of processs allowed to be used since "wait" was just printed
        $process_batches_count = $process_batches_count + 1;
    }
    return $process_batches_count;
}

1;
