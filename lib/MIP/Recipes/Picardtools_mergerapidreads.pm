package MIP::Recipes::Picardtools_mergerapidreads;

use strict;
use warnings;
use warnings qw{FATAL utf8};
use utf8;
use open qw{:encoding(UTF-8) :std};
use autodie qw{:all};
use charnames qw{:full :short};
use Carp;
use English qw{-no_match_vars};
use Params::Check qw{check allow last_error};
use File::Basename qw(dirname fileparse);
use File::Spec::Functions qw(catdir catfile devnull);

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{analysis_bwa_mem};

}

sub picardtools_mergerapidreads {

##picardtools_mergerapidreads

##Function : Merges all batch read processes to one file using Picardtools mergesamfiles within each sampleid. The read batch proccessed files have to be sorted before attempting to merge.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id, $program_name, $outaligner_dir_ref, $temp_directory_ref
##         : $parameter_href             => Parameter hash {REF}
##         : $active_parameter_href      => Active parameters for this analysis hash {REF}
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $file_info_href             => The file_info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href                => Job id hash {REF}
##         : $sample_id_ref              => Sample id {REF}
##         : $program_name               => Program name
##         : $outaligner_dir_ref         => Outaligner_dir used in the analysis {REF}
##         : $temp_directory_ref         => Temporary directory {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $temp_directory_ref;
    my $outaligner_dir_ref;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id_ref;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
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
        sample_id_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$sample_id_ref
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw(setup_script);
    use MIP::Processmanagement::Processes qw(print_wait);
    use MIP::Language::Java qw{java_core};
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $$sample_id_ref,
            program_name          => $program_name,
            program_directory     => lc($$outaligner_dir_ref),
            core_number  => $active_parameter_href->{max_cores_per_node},
            process_time => 20,
        }
    );

    ## Assign directories
    my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
        $$sample_id_ref, $$outaligner_dir_ref );
    my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
        $$sample_id_ref, $$outaligner_dir_ref );

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$$sample_id_ref}{pbwa_mem}{file_tag};
    my $outfile_tag =
      $file_info_href->{$$sample_id_ref}{ "p" . $program_name }{file_tag};

    my $process_batches_count = 1;
    my $core_tracker          = 0
      ; #Required to portion out cores and files before wait and to track the MOS_BU outfiles to correct lane

    for (
        my $infile_counter = 0 ;
        $infile_counter <
        scalar( @{ $infile_lane_prefix_href->{$$sample_id_ref} } ) ;
        $infile_counter++
      )
    {    #For all files from

        my $infile =
          $infile_lane_prefix_href->{$$sample_id_ref}[$infile_counter];
        my $nr_read_batch_process =
          $sample_info_href->{sample}{$$sample_id_ref}
          { $infile_lane_prefix_href->{$$sample_id_ref}[$infile_counter] }
          {pbwa_mem}{read_batch_process};

        if ( $nr_read_batch_process > 0 )
        {    #Check that we have read batch processes to merge

            $process_batches_count = print_wait(
                {
                    process_counter => $core_tracker,
                    max_process_number =>
                      $active_parameter_href->{max_cores_per_node},
                    process_batches_count => $process_batches_count,
                    FILEHANDLE            => $FILEHANDLE,
                }
            );

            for (
                my $read_batch_processes_count = 0 ;
                $read_batch_processes_count < $nr_read_batch_process ;
                $read_batch_processes_count++
              )
            {

                if ( $read_batch_processes_count eq 0 ) {

                    java_core(
                        {
                            FILEHANDLE        => $FILEHANDLE,
                            memory_allocation => "Xmx4g",
                            java_use_large_pages =>
                              $active_parameter_href->{java_use_large_pages},
                            temp_directory =>
                              $active_parameter_href->{temp_directory},
                            java_jar => catfile(
                                $active_parameter_href->{picardtools_path},
                                "picard.jar"
                            ),
                        }
                    );

                    print $FILEHANDLE "MergeSamFiles ";
                    print $FILEHANDLE "USE_THREADING=TRUE "
                      ; #Create a background thread to encode, compress and write to disk the output file
                    print $FILEHANDLE "CREATE_INDEX=TRUE "
                      ; #Create a BAM index when writing a coordinate-sorted BAM file.
                    print $FILEHANDLE "OUTPUT="
                      . catfile( $outsample_directory,
                        $infile_lane_prefix_href->{$$sample_id_ref}
                          [$infile_counter] . $outfile_tag . q{.bam} )
                      . " ";    #OutFile
                }
                print $FILEHANDLE "INPUT="
                  . catfile( $insample_directory,
                    $infile_lane_prefix_href->{$$sample_id_ref}[$infile_counter]
                      . "_"
                      . $read_batch_processes_count
                      . $outfile_tag
                      . q{.bam} )
                  . " ";        #InFile(s)
            }
            say $FILEHANDLE "& ", "\n";
            $core_tracker++
              ; #Track nr of merge calls for infiles so that wait can be printed at the correct intervals (dependent on $active_parameter_href->{max_cores_per_node})
        }
        else
        { #Still needs to rename file to be included in potential merge of BAM files in next step

            java_core(
                {
                    FILEHANDLE        => $FILEHANDLE,
                    memory_allocation => "Xmx4g",
                    java_use_large_pages =>
                      $active_parameter_href->{java_use_large_pages},
                    temp_directory => $active_parameter_href->{temp_directory},
                    java_jar       => catfile(
                        $active_parameter_href->{picardtools_path},
                        "picard.jar"
                    ),
                }
            );

            print $FILEHANDLE "MergeSamFiles ";
            print $FILEHANDLE "USE_THREADING=TRUE "
              ; #Create a background thread to encode, compress and write to disk the output file
            print $FILEHANDLE "CREATE_INDEX=TRUE "
              ;   #Create a BAM index when writing a coordinate-sorted BAM file.
            print $FILEHANDLE "INPUT="
              . catfile( $insample_directory,
                    $infile_lane_prefix_href->{$$sample_id_ref}[$infile_counter]
                  . "_0"
                  . $outfile_tag
                  . "_rg.bam" )
              . " ";    #InFile
            say $FILEHANDLE "OUTPUT="
              . catfile( $outsample_directory,
                    $infile_lane_prefix_href->{$$sample_id_ref}[$infile_counter]
                  . $outfile_tag
                  . q{.bam} )
              . " &";    #OutFile
        }
    }
    say $FILEHANDLE q{wait}, "\n";

    ## Remove temp directory
    gnu_rm(
        {
            infile_path => $active_parameter_href->{temp_directory},
            force       => 1,
            recursive   => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );

    close($FILEHANDLE);

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                job_id_href             => $job_id_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_ids_ref => \@{ $active_parameter_href->{sample_ids} },
                family_id      => $active_parameter_href->{family_id},
                path => $parameter_href->{ "p" . $program_name }{chain},
                log  => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
}

1;
