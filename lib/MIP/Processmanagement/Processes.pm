package MIP::Processmanagement::Processes;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    # Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{check allow last_error};

use FindBin qw{$Bin};    # Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

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
    our @EXPORT_OK =
      qw{add_to_job_id_dependency_string add_parallel_job_id_to_dependency_tree add_pan_job_id_to_sample_id_dependency_tree add_pan_job_id_to_family_id_dependency_tree add_sample_job_id_to_dependency_tree print_wait};

}

##Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};
Readonly my $EMPTY_STR  => q{};
Readonly my $COLON      => q{:};

sub add_to_job_id_dependency_string {

##add_to_job_id_dependency_string

##Function : Adds all previous jobIds per familyChainKey and chain_key to job_ids dependency string used to set the dependency in SLURM.
##Returns  : "$job_ids"
##Arguments: $job_id_href, $family_id_chain_key, $chain_key
##         : $job_id_href         => The info on jobIds hash {REF}
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

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    # Job_id string to submit to workload manager
    my $job_ids_string =
      $EMPTY_STR;

    # Alias hash to get job_id dependencies from
    my $chain_job_id_href = $job_id_href->{$family_id_chain_key}{$chain_key};

    if ($chain_job_id_href) {

      JOB_IDS:
        while ( my ( $job_index, $job_id ) = each @{$chain_job_id_href} ) {

	  # Only for defined job_ids
	  if (defined $job_id) {

            # Only 1 previous job_id
            if ( ( !$job_index ) && ( scalar @{$chain_job_id_href} == 1 ) ) {

                #Single job_id start with $COLON and end without $COLON
                $job_ids_string .= $COLON . $job_id;
            }
            elsif ( !$job_index ) {
                #First job_id to add

                #First job_id start with $COLON and ends with $COLON
                $job_ids_string .= $COLON . $job_id . $COLON;
            }
            elsif ( $job_index eq ( scalar @{$chain_job_id_href} - 1 ) ) {
                #Last job_id

                #Last job_id finish without $COLON
                $job_ids_string .= $job_id;
            }
            else {
                #JobIDs in the middle

                $job_ids_string .= $job_id . $COLON;
            }
	  }
        }
    }
    return $job_ids_string;
}

sub add_parallel_job_id_to_dependency_tree {

##add_parallel_job_id_to_dependency_tree

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

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $chain_key_type = q{parallel};
    my $parallel_jobs_chain_key;

    ## Push parallel job_ids
  INFILES:
    while ( my ($infile_index) = each $infile_lane_prefix_href->{$sample_id} ) {

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

            push $sample_id_job_ids_ref, @{$job_ids_ref};
        }
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

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $pan_chain_key;

    ## Push pan jobs

    # Set key
    $pan_chain_key =
      $family_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    # If pan job ids exists
    if ( exists $job_id_href->{$family_id_chain_key}{$pan_chain_key} ) {

        # Alias pan job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}{$pan_chain_key} };

        # Alias job_ids array for sample_id in job_id_href to push to
        my $sample_id_job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} };

        # Add job_ids to sample_id job_ids
        push $sample_id_job_ids_ref, @{$job_ids_ref};
    }
    return;
}

sub add_pan_job_id_to_family_id_dependency_tree {

##add_pan_job_id_to_family_id_dependency_tree

##Function : Saves pan (i.e job_ids that affect all chains) job_id to the the family_id chain.
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

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $pan_chain_key;

    ## Push family_pan jobs

    # Set key
    $pan_chain_key = $family_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

    # If family_pan job ids exists
    if ( exists $job_id_href->{$family_id_chain_key}{$pan_chain_key} ) {

        # Alias family_pan job_ids array to push from
        my $job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}{$pan_chain_key} };

        ## Use $family_id_chain_key instead of $sample_id_chain_key
        # Alias job_ids array for family_id in job_id_href to push to
        my $family_id_job_ids_ref =
          \@{ $job_id_href->{$family_id_chain_key}{$family_id_chain_key} };

        # Add job_ids to family_id job_ids
        push $family_id_job_ids_ref, @{$job_ids_ref};
    }
    return;
}

sub add_sample_job_id_to_dependency_tree {

##add_sample_job_id_to_dependency_tree

##Function : Saves job_id to the correct hash array depending on chain type.
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

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Push sample job

    # Alias job_ids array for sample_id in job_id_href to push to
    my $sample_id_job_ids_ref =
      \@{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key} };

    # Add to sample_id job dependency tree
    push $sample_id_job_ids_ref, $job_id_returned;

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

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

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
