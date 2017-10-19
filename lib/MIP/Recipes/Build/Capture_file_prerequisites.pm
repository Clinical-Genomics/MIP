package MIP::Recipes::Build::Capture_file_prerequisites;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_capture_file_prerequisites };

}

##Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub build_capture_file_prerequisites {

## Function : Creates the target "infiles_list" "padded.infile_list" and interval_list files.
## Returns  :
## Arguments: $parameter_href              => Parameter hash {REF}
##          : $active_parameter_href       => Active parameters for this analysis hash {REF}
##          : $sample_info_href            => Info on samples and family hash {REF}
##          : $infile_lane_prefix_href     => Infile(s) without the ".ending" {REF}
##          : $job_id_href                 => Job id hash {REF}
##          : $infile_list_suffix          => Infile list suffix
##          : $padded_infile_list_suffix   => Padded infile list suffix
##          : $padded_interval_list_suffix => Padded interval list suffix
##          : $program_name                => Program name
##          : $FILEHANDLE                  => Filehandle to write to
##          : $log                         => Log object
##          : $family_id                   => Family ID
##          : $outaligner_dir              => Outaligner_dir used in the analysis {REF}
##          : $temp_directory              => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $infile_list_suffix;
    my $padded_infile_list_suffix;
    my $padded_interval_list_suffix;
    my $program_name;
    my $FILEHANDLE;
    my $log;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;

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
        infile_list_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_list_suffix,
        },
        padded_infile_list_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$padded_infile_list_suffix,
        },
        padded_interval_list_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$padded_interval_list_suffix,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        FILEHANDLE => { store => \$FILEHANDLE, },
        log        => { store => \$log, },
        family_id  => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Gnu::Coreutils qw{ gnu_rm gnu_cat gnu_ln };
    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Program::Fasta::Picardtools
      qw{ picardtools_createsequencedictionary };
    use MIP::Program::Interval::Picardtools qw{ picardtools_intervallisttools };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_add_to_samples };

    ## Constants
    Readonly my $NR_OF_BASES_PADDING => 100;
    Readonly my $MAX_RANDOM_NUMBER   => 100_00;

    my $file_path;

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

    ## No supplied FILEHANDLE i.e. create new sbatch script
    if ( not defined $FILEHANDLE ) {

        ## Create anonymous filehandle
        $FILEHANDLE = IO::Handle->new();

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

    ## Generate a random integer.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

  BED_FILE:
    foreach my $exome_target_bed_file (
        keys %{ $active_parameter_href->{exome_target_bed} } )
    {

        $log->warn( q{Will try to create required }
              . $exome_target_bed_file
              . q{ associated file(s) before executing }
              . $program_name );

        ## Add random integer
        my $exome_target_bed_file_random =
          $exome_target_bed_file . $UNDERSCORE . $random_integer;

        say {$FILEHANDLE} q{## CreateSequenceDictionary from reference};

        picardtools_createsequencedictionary(
            {
                FILEHANDLE        => $FILEHANDLE,
                memory_allocation => q{Xmx2g},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $temp_directory,
                java_jar       => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
                referencefile_path => $referencefile_path,
                outfile_path => $exome_target_bed_file_random . $DOT . q{dict},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        say {$FILEHANDLE}
          q{## Add target file to headers from sequence dictionary};
        gnu_cat(
            {
                infile_paths_ref => [
                    $exome_target_bed_file_random . $DOT . q{dict},
                    $exome_target_bed_file
                ],
                outfile_path => $exome_target_bed_file_random
                  . $DOT
                  . q{dict_body},
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Remove track and browser header lines and reformat columns in file
        _reformat_capture_file(
            {
                exome_target_bed_file_random => $exome_target_bed_file_random,
                FILEHANDLE                   => $FILEHANDLE,
            }
        );

        say {$FILEHANDLE} q{## Create} . $infile_list_suffix;

        my @infile_paths_ref =
          (     $exome_target_bed_file_random
              . $DOT
              . q{dict_body_col_5.interval_list} );
        my $infile_list_outfile_path =
            $exome_target_bed_file_random
          . $DOT
          . q{dict_body_col_5}
          . $UNDERSCORE
          . $infile_list_suffix;
        picardtools_intervallisttools(
            {
                memory_allocation => q{Xmx2g},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $temp_directory,
                java_jar       => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
                infile_paths_ref   => \@infile_paths_ref,
                outfile_path       => $infile_list_outfile_path,
                FILEHANDLE         => $FILEHANDLE,
                referencefile_path => $referencefile_path,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        my $intended_file_path = $exome_target_bed_file . $infile_list_suffix;
        my $temporary_file_path =
            $exome_target_bed_file_random
          . $DOT
          . q{dict_body_col_5}
          . $UNDERSCORE
          . $infile_list_suffix;

        ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
        check_exist_and_move_file(
            {
                FILEHANDLE          => $FILEHANDLE,
                intended_file_path  => $intended_file_path,
                temporary_file_path => $temporary_file_path,
            }
        );

        say {$FILEHANDLE} q{#Create} . $padded_infile_list_suffix;

        my $padded_infile_list_outfile_path =
            $exome_target_bed_file_random
          . $DOT
          . q{dict_body_col_5}
          . $padded_infile_list_suffix;
        picardtools_intervallisttools(
            {
                memory_allocation => q{Xmx2g},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $active_parameter_href->{temp_directory},
                java_jar       => catfile(
                    $active_parameter_href->{picardtools_path},
                    q{picard.jar}
                ),
                infile_paths_ref   => \@infile_paths_ref,
                outfile_path       => $padded_infile_list_outfile_path,
                padding            => $NR_OF_BASES_PADDING,
                FILEHANDLE         => $FILEHANDLE,
                referencefile_path => $referencefile_path,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        $intended_file_path =
          $exome_target_bed_file . $padded_infile_list_suffix;
        $temporary_file_path =
            $exome_target_bed_file_random
          . $DOT
          . q{dict_body_col_5}
          . $padded_infile_list_suffix;

        ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
        check_exist_and_move_file(
            {
                FILEHANDLE          => $FILEHANDLE,
                intended_file_path  => $intended_file_path,
                temporary_file_path => $temporary_file_path,
            }
        );

        say {$FILEHANDLE} q{#Create }
          . $padded_interval_list_suffix
          . q{ by softlinking};

        gnu_ln(
            {
                force       => 1,
                symbolic    => 1,
                target_path => $exome_target_bed_file
                  . $padded_infile_list_suffix,
                link_path => $exome_target_bed_file
                  . $padded_interval_list_suffix,
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Remove temporary files
        say {$FILEHANDLE} q{#Remove temporary files};

        my @temp_files = (
            $exome_target_bed_file_random
              . $DOT
              . q{dict_body_col_5.interval_list},
            $exome_target_bed_file_random . $DOT . q{dict_body},
            $exome_target_bed_file_random . $DOT . q{dict},
        );
      FILE_TO_REMOVE:
        foreach my $file (@temp_files) {

            gnu_rm(
                {
                    infile_path => $file,
                    force       => 1,
                    FILEHANDLE  => $FILEHANDLE,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }
    ## Unless FILEHANDLE was supplied close filehandle and submit
    if ( not defined $FILEHANDLE ) {

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

sub _reformat_capture_file {

## Function : Remove track and browser header lines and reformat columns in file
## Returns  :
## Arguments: $exome_target_bed_file_random => Exome target bed file
##          : $FILEHANDLE                   => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_file_random;
    my $FILEHANDLE;

    my $tmpl = {
        exome_target_bed_file_random => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$exome_target_bed_file_random,
        },
        FILEHANDLE => { store => \$FILEHANDLE, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Remove unnecessary info and reformat
    say {$FILEHANDLE}
      q{#Remove target annotations, 'track', 'browse' and keep only 5 columns};

    ## Execute perl
    print {$FILEHANDLE} q?perl -nae '?;

    ## If header line
    print {$FILEHANDLE} q?if ($_=~/@/) { ?;

    ## Write to stdout
    print {$FILEHANDLE} q?print $_;} ?;

    ## If track - do nothing
    print {$FILEHANDLE} q?elsif ($_=~/^track/) {} ?;

    ## If browser - do nothing
    print {$FILEHANDLE} q?elsif ($_=~/^browser/) {} ?;

    ## Else print reformated line to stdout
    print {$FILEHANDLE}
q?else {print @F[0], "\t", (@F[1] + 1), "\t", @F[2], "\t", "+", "\t", "-", "\n";}' ?;

    ## Infile
    print {$FILEHANDLE} $exome_target_bed_file_random
      . $DOT
      . q{dict_body}
      . $SPACE;

    ## Write to
    print {$FILEHANDLE} q{>} . $SPACE;
    say   {$FILEHANDLE} $exome_target_bed_file_random . $DOT
      . q{dict_body_col_5.interval_list}, $NEWLINE;

    return;
}

1;
