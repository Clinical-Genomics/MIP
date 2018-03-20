package MIP::Recipes::Build::Human_genome_prerequisites;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_human_genome_prerequisites };

}

##Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub build_human_genome_prerequisites {

## Function : Creates the human genome prerequisites using active_parameters{human_genome_reference} as reference.
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $program_name            => Program under evaluation
##          : $FILEHANDLE              => Filehandle to write to. A new sbatch script will be generated if $FILEHANDLE is lacking, else write to exising $FILEHANDLE {Optional}
##          : $log                     => Log object
##          : $random_integer          => The random integer to create temporary file name
##          : $family_id               => Family ID
##          : $reference_dir           => MIP reference directory
##          : $outaligner_dir          => The outaligner_dir used in the analysis
##          : $human_genome_reference  => Human genome reference

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $FILEHANDLE;
    my $log;
    my $random_integer;

    ## Default(s)
    my $family_id;
    my $reference_dir;
    my $outaligner_dir;
    my $human_genome_reference;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        FILEHANDLE     => { store       => \$FILEHANDLE, },
        log            => { store       => \$log },
        random_integer => { strict_type => 1, store => \$random_integer, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        human_genome_reference => {
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Reference qw{ check_capture_file_prerequisites };
    use MIP::Gnu::Coreutils qw{ gnu_rm gnu_ln };
    use MIP::Language::Java qw{ java_core };
    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Program::Alignment::Samtools qw{ samtools_faidx };
    use MIP::Program::Compression::Gzip qw{ gzip };
    use MIP::Program::Fasta::Picardtools
      qw{ picardtools_createsequencedictionary };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_add_to_samples };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 100_00;

    my $file_path;
    my $submit_switch;

    ## Alias
    my $mip_program_mode = $active_parameter_href->{ q{p} . $program_name };

    ## No supplied FILEHANDLE i.e. create new sbatch script
    if ( not defined $FILEHANDLE ) {

        $submit_switch = 1;

        ## Create anonymous filehandle
        $FILEHANDLE = IO::Handle->new();

        ## Generate a random integer between 0-10,000.
        $random_integer = int rand $MAX_RANDOM_NUMBER;

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ($file_path) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $family_id,
                program_name          => $program_name,
                program_directory     => $outaligner_dir,
            }
        );
    }

    ## Check for compressed files
    if ( $file_info_href->{human_genome_compressed} eq q{compressed} ) {

        $log->warn( q{Will try to decompress }
              . $human_genome_reference
              . q{ before executing }
              . $program_name );

        ## Perl wrapper for writing gzip recipe to $FILEHANDLE
        gzip(
            {
                decompress  => 1,
                infile_path => $human_genome_reference,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Replace the .fasta.gz ending with .fasta since this will execute before the analysis
        ## Hence, changing the original file name ending from ".fastq" to ".fastq.gz".
        $human_genome_reference =~ s/.fasta.gz/.fasta/xsmg;

        $log->info(
            q{Set human_genome_reference to: } . $human_genome_reference );
        $file_info_href->{human_genome_compressed} = q{uncompressed};
    }

    if ( exists $file_info_href->{exome_target_bed} ) {

        check_capture_file_prerequisites(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                infile_list_suffix => $file_info_href->{exome_target_bed}[0],
                padded_infile_list_suffix =>
                  $file_info_href->{exome_target_bed}[1],
                padded_interval_list_suffix =>
                  $file_info_href->{exome_target_bed}[2],
                program_name => $program_name,
                FILEHANDLE   => $FILEHANDLE,
                log          => $log,
            }
        );
    }
    if ( $parameter_href->{q{human_genome_reference}}{build_file} == 1 ) {

      FILE_ENDING:
        foreach my $file_ending (
            @{ $file_info_href->{human_genome_reference_file_endings} } )
        {

            if ( $file_ending eq $DOT . q{dict} ) {

                $log->warn( q{Will try to create }
                      . $file_ending
                      . q{ file for }
                      . $human_genome_reference
                      . q{ before executing }
                      . $program_name );

                my $filename_prefix = catfile( $reference_dir,
                    $file_info_href->{human_genome_reference_name_prefix} );

                say {$FILEHANDLE} q{#CreateSequenceDictionary from reference};

                picardtools_createsequencedictionary(
                    {
                        FILEHANDLE        => $FILEHANDLE,
                        memory_allocation => q{Xmx2g},
                        java_use_large_pages =>
                          $active_parameter_href->{java_use_large_pages},
                        temp_directory =>
                          $active_parameter_href->{temp_directory},
                        java_jar => catfile(
                            $active_parameter_href->{picardtools_path},
                            q{picard.jar}
                        ),
                        referencefile_path => $human_genome_reference,
                        outfile_path       => $filename_prefix
                          . $UNDERSCORE
                          . $random_integer
                          . $file_ending,
                    }
                );
                say {$FILEHANDLE} $NEWLINE;

                my $intended_file_path = $filename_prefix . $file_ending;
                my $temporary_file_path =
                    $filename_prefix
                  . $UNDERSCORE
                  . $random_integer
                  . $file_ending;

                ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
                check_exist_and_move_file(
                    {
                        FILEHANDLE          => $FILEHANDLE,
                        intended_file_path  => $intended_file_path,
                        temporary_file_path => $temporary_file_path,
                    }
                );
            }
            if ( $file_ending eq $DOT . q{fai} ) {

                $log->warn( q{Will try to create }
                      . $file_ending
                      . q{ file for }
                      . $human_genome_reference
                      . q{ before executing }
                      . $program_name );

                my $human_genome_reference_temp_file =
                  $human_genome_reference . $UNDERSCORE . $random_integer;

                say {$FILEHANDLE} q{## Fai file from reference};
                gnu_ln(
                    {
                        force       => 1,
                        symbolic    => 1,
                        target_path => $human_genome_reference,
                        link_path   => $human_genome_reference_temp_file,
                        FILEHANDLE  => $FILEHANDLE,
                    }
                );
                say {$FILEHANDLE} $NEWLINE;

                samtools_faidx(
                    {
                        infile_path => $human_genome_reference_temp_file,
                        FILEHANDLE  => $FILEHANDLE,
                    }
                );
                say {$FILEHANDLE} $NEWLINE;

                my $intended_file_path = $human_genome_reference . $file_ending;
                my $temporary_file_path =
                  $human_genome_reference_temp_file . $file_ending;

                ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
                check_exist_and_move_file(
                    {
                        FILEHANDLE          => $FILEHANDLE,
                        intended_file_path  => $intended_file_path,
                        temporary_file_path => $temporary_file_path,
                    }
                );

                ## Remove soft link
                gnu_rm(
                    {
                        infile_path => $human_genome_reference_temp_file,
                        force       => 1,
                        FILEHANDLE  => $FILEHANDLE,
                    }
                );
                say {$FILEHANDLE} $NEWLINE;
            }
        }

        ## Only create once
        $parameter_href->{human_genome_reference}{build_file} = 0;
    }

    ## Unless FILEHANDLE was supplied close it and submit
    if ($submit_switch) {

        close $FILEHANDLE;

        if ( $mip_program_mode == 1 ) {

            slurm_submit_job_no_dependency_add_to_samples(
                {
                    job_id_href => $job_id_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $family_id,
                    path             => q{MAIN},
                    sbatch_file_name => $file_path,
                    log              => $log,
                }
            );
        }
    }
    return;
}

1;
