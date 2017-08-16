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
    our @EXPORT_OK = qw{push_to_job_id print_wait};

}

##Constants
Readonly my $NEWLINE => qq{\n};

sub push_to_job_id {

##push_to_job_id

##Function : Saves job_id to the correct hash array depending on chain type.
##Returns  : ""
##Arguments: $job_id_href, $infile_lane_prefix_href, $family_id_chain_key, $sample_id_chain_key, $sample_id, $path, $chain_key_type
##         : $job_id_href             => Info on jobIds hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $family_id_chain_key     => Family ID chain hash key
##         : $sample_id_chain_key     => Sample ID chain hash key
##         : $sample_id               => Sample ID
##         : $path                    => Trunk or branch
##         : $chain_key_type          => 'parallel', 'merged' or 'family_merged'

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parallel_chains_ref;
    my $family_id_chain_key;
    my $sample_id_chain_key;
    my $sample_id;
    my $path;
    my $chain_key_type;

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
        parallel_chains_ref =>
          { default => [], strict_type => 1, store => \$parallel_chains_ref },
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
        chain_key_type => {
            required    => 1,
            defined     => 1,
            allow       => [qw{parallel merged family_merged}],
            strict_type => 1,
            store       => \$chain_key_type
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $chain_key;

    ## Push parallel job_ids
    if ( $chain_key_type eq q{parallel} ) {

      INFILES:
        while ( my ($infile_index) =
            each $infile_lane_prefix_href->{$sample_id} )
        {

            # Set key
            $chain_key =
                $sample_id . q{_}
              . $chain_key_type . q{_}
              . $path
              . $infile_index;

            ## If parellel job_ids exists
            if ( exists $job_id_href->{$family_id_chain_key}{$chain_key} ) {

                # Alias job_ids array for sample_id in job_id_href to push to
                my $sample_id_job_ids_ref =
                  \@{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key}
                  };

                # Alias parallel job_ids array to push from
                my $job_ids_ref =
                  \@{ $job_id_href->{$family_id_chain_key}{$chain_key} };

                push $sample_id_job_ids_ref, @{$job_ids_ref};
            }
        }
    }
    elsif (( $chain_key_type eq q{merged} )
        || ( $chain_key_type eq q{family_merged} ) )
    {
        # Push merged jobs

        # Set key
        $chain_key = $family_id_chain_key . q{_} . $sample_id_chain_key;

        # If merged job ids exists
        if ( exists $job_id_href->{$family_id_chain_key}{$chain_key} ) {

            # Alias parallel job_ids array to push from
            my $job_ids_ref =
              \@{ $job_id_href->{$family_id_chain_key}{$chain_key} };

            # Use $family_id_chain_key instead of $sample_id_chain_key
            if ( $chain_key_type eq q{family_merged} ) {

                # Alias job_ids array for sample_id in job_id_href to push to
                my $family_id_job_ids_ref =
                  \@{ $job_id_href->{$family_id_chain_key}{$family_id_chain_key}
                  };

                # Add job_ids to family_id job_ids
                push $family_id_job_ids_ref, @{$job_ids_ref};
            }
            else {

                # Alias job_ids array for sample_id in job_id_href to push to
                my $sample_id_job_ids_ref =
                  \@{ $job_id_href->{$family_id_chain_key}{$sample_id_chain_key}
                  };

                # Add job_ids to sample_id job_ids
                push $sample_id_job_ids_ref, @{$job_ids_ref};
            }
        }
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
