package MIP::Recipes::Markduplicates;

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
use File::Spec::Functions qw{ catdir catfile devnull };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_markduplicates };

}

##Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_markduplicates {

## analysis_markduplicates

## Function : Mark duplicated reads using Picardtools markduplicates or Sambamba markduplicates in files generated from alignment (sorted, merged).
## Returns  : |$xargs_file_counter
## Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $sample_id, $program_name, $program_info_path, $file_path, $FILEHANDLE, family_id, $temp_directory, $outaligner_dir, $xargs_file_counter
##          : $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $sample_id           => Sample id {REF}
##          : $program_name            => Program name
##          : $program_info_path       => The program info path
##          : $file_path               => File path
##          : $FILEHANDLE              => Filehandle to write to
##          : $family_id               => Family id
##          : $temp_directory          => Temporary directory
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $sample_id;
    my $program_name;
    my $program_info_path;
    my $file_path;
    my $FILEHANDLE;

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
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        FILEHANDLE    => { store => \$FILEHANDLE },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        outaligner_dir => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Recipes::Xargs qw{ xargs_command };
    use MIP::Program::Alignment::Sambamba qw{ sambamba_flagstat };
    use MIP::Program::Alignment::Picardtools qw{ picardtools_markduplicates };
    use MIP::Gnu::Coreutils qw{ gnu_cat };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_sample };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};
    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    unless ( defined $FILEHANDLE ) {           #Run as individual sbatch script

        $FILEHANDLE = IO::Handle->new();        #Create anonymous filehandle
    }

    ## Assign directories
    my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
        $sample_id, $outaligner_dir );
    my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
        $sample_id, $outaligner_dir );

    # Used downstream
    $parameter_href->{$mip_program_name}{$sample_id}{indirectory} =
      $outsample_directory;

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{ppicardtools_mergesamfiles}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$mip_program_name}{file_tag};

    ## Files
    my $infile_prefix       = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix      = $merged_infile_prefix . $outfile_tag;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
    my $infile_suffix = my $outfile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
            jobid_chain    => $job_id_chain,
        }
    );

    ## Sums all mapped and duplicate reads and takes fraction of before finishing
    my $regexp =
q?perl -nae'my %feature; while (<>) { if($_=~/duplicates/ && $_=~/^(\d+)/) {$feature{dup} = $feature{dup} + $1} if($_=~/\d+\smapped/ && $_=~/^(\d+)/) {$feature{map} = $feature{map} + $1} } print "Read Mapped: ".$feature{map}."\nDuplicates: ".$feature{dup}."\n"."Fraction Duplicates: ".$feature{dup}/$feature{map}, "\n"; last;'?;

    # Store which program performed the markduplication
    my $markduplicates_program;

    if ( not $$reduce_io_ref ) {    #Run as individual sbatch script

        use MIP::Script::Setup_script qw{ setup_script };

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $sample_id,
                program_name          => $program_name,
                program_directory     => lc($outaligner_dir),
                core_number           => $core_number,
                process_time => $time,
                temp_directory => $temp_directory
            }
        );

        ## Copy file(s) to temporary directory
        say $FILEHANDLE q{## Copy file(s) to temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE      => $FILEHANDLE,
                XARGSFILEHANDLE => $XARGSFILEHANDLE,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                file_path   => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                infile             => $infile_prefix,
                indirectory        => $insample_directory,
                file_ending        => substr( $infile_suffix, 0, 2 )
                  . "*",    #q{.bam} -> ".b*" for getting index as well
                temp_directory => $temp_directory,
            }
        );
    }

    ## Marking Duplicates
    say $FILEHANDLE q{## Marking Duplicates};

    ##Picardtools
    if ( $active_parameter_href->{markduplicates_picardtools_markduplicates} ) {

        $markduplicates_program = q{picardtools_markduplicates};

        ## Create file commands for xargs
        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                first_command      => q{java},
                memory_allocation  => q{Xmx4g},
                java_use_large_pages =>
                  $active_parameter_href->{java_use_large_pages},
                temp_directory => $temp_directory,
                java_jar       => catfile( $active_parameter_href->{picardtools_path}
                  , q{picard.jar}),
            }
        );

      CONTIG:
        foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

            picardtools_markduplicates(
                {
                    infile_paths_ref =>
                      [ $file_path_prefix . "_" . $contig . $infile_suffix ],
                    outfile_path => $outfile_path_prefix . "_"
                      . $contig
                      . $outfile_suffix,
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $contig
                      . ".stderr.txt",
                    metrics_file => $outfile_path_prefix . "_"
                      . $contig
                      . ".metric",
                    FILEHANDLE => $XARGSFILEHANDLE,
                    referencefile_path =>
                      $active_parameter_href->{human_genome_reference},
                    create_index => "true",
                }
            );
            print {$XARGSFILEHANDLE} q{;} . $SPACE;

            ## Process BAM with sambamba flagstat to produce metric file for downstream analysis
            sambamba_flagstat(
                {
                    infile_path => $outfile_path_prefix . "_"
                      . $contig
                      . $outfile_suffix,
                    outfile_path => $outfile_path_prefix . "_"
                      . $contig
                      . "_metric",
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $contig
                      . ".stderr.txt",
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );
            say $XARGSFILEHANDLE $NEWLINE;
        }
    }

    ## Sambamba
    if ( $active_parameter_href->{markduplicates_sambamba_markdup} ) {

        $markduplicates_program = q{sambamba_markdup};

        use MIP::Program::Alignment::Sambamba qw(sambamba_markdup);

        ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
            {
                FILEHANDLE         => $FILEHANDLE,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                file_path          => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
            }
        );

      CONTIG:
        foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

            sambamba_markdup(
                {
                    infile_path => $file_path_prefix . "_"
                      . $contig
                      . $infile_suffix,
                    outfile_path => $outfile_path_prefix . "_"
                      . $contig
                      . $outfile_suffix,
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $contig
                      . ".stderr.txt",
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    temp_directory  => $temp_directory,
                    hash_table_size => $active_parameter_href
                      ->{markduplicates_sambamba_markdup_hash_table_size},
                    overflow_list_size => $active_parameter_href
                      ->{markduplicates_sambamba_markdup_overflow_list_size},
                    io_buffer_size => $active_parameter_href
                      ->{markduplicates_sambamba_markdup_io_buffer_size},
                    show_progress => 1,
                }
            );
            print {$XARGSFILEHANDLE} q{;} . $SPACE;

            ## Process BAM with sambamba flagstat to produce metric file for downstream analysis
            sambamba_flagstat(
                {
                    infile_path => $outfile_path_prefix . "_"
                      . $contig
                      . $outfile_suffix,
                    outfile_path => $outfile_path_prefix . "_"
                      . $contig
                      . "_metric",
                    stderrfile_path => $xargs_file_path_prefix . "."
                      . $contig
                      . ".stderr.txt",
                    FILEHANDLE => $XARGSFILEHANDLE,
                }
            );
            say $XARGSFILEHANDLE $NEWLINE;
        }
    }

    ## Concatenate all metric files
    gnu_cat(
        {
            infile_paths_ref => [ $outfile_path_prefix . "_*_metric" ],
            outfile_path     => $outfile_path_prefix . "_metric_all",
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    say $FILEHANDLE $NEWLINE;

    ## Sum metric over concatenated file
    print $FILEHANDLE $regexp . $SPACE;
    print $FILEHANDLE $outfile_path_prefix . "_metric_all" . $SPACE;
    say $FILEHANDLE "> " . $outfile_path_prefix . "_metric" . $SPACE,
      $NEWLINE;

    migrate_file(
        {
            infile_path  => $outfile_path_prefix . q{_metric},
            outfile_path => $outsample_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say $FILEHANDLE q{wait}, $NEWLINE;

    if ( $mip_program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                sample_id        => $sample_id,
                program_name     => q{markduplicates},
                infile           => $merged_infile_prefix,
                outdirectory     => $outsample_directory,
                outfile          => $outfile_prefix . q{_metric},
            }
        );

# Markduplicates can be processed by either picardtools markduplicates or sambamba markdup
        $sample_info_href->{sample}{$sample_id}{program}{markduplicates}
          {$merged_infile_prefix}{processed_by} = $markduplicates_program;

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            my $most_complete_format_key =
              q{most_complete_} . substr( $outfile_suffix, 1 );
            $sample_info_href->{sample}{$sample_id}
              {$most_complete_format_key}{path} = catfile(
                $outsample_directory,
                $outfile_prefix . "_"
                  . $file_info_href->{contigs_size_ordered}[0]
                  . $outfile_suffix
              );
        }
    }

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        ## Copies file from temporary directory. Per contig
        say $FILEHANDLE q{## Copy file from temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                FILEHANDLE      => $FILEHANDLE,
                XARGSFILEHANDLE => $XARGSFILEHANDLE,
                contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
                file_path   => $file_path,
                program_info_path  => $program_info_path,
                core_number        => $core_number,
                xargs_file_counter => $xargs_file_counter,
                outfile            => $outfile_prefix,
                outdirectory       => $outsample_directory,
                temp_directory     => $temp_directory,
                file_ending        => substr( $infile_suffix, 0, 2 ) . "*",
            }
        );
    }
    else {

        ## Remove file at temporary Directory
        remove_contig_files(
            {
                file_elements_ref =>
                  \@{ $file_info_href->{contigs_size_ordered} },
                FILEHANDLE  => $FILEHANDLE,
                core_number => $core_number,
                file_name   => $infile_prefix,
                file_ending => substr( $infile_suffix, 0, 2 ) . "*",
                indirectory => $temp_directory,
            }
        );
    }

    close $XARGSFILEHANDLE;

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        close $FILEHANDLE;

        if ( $mip_program_mode == 1 ) {

            slurm_submit_job_sample_id_dependency_add_to_sample(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    family_id               => $family_id,
                    sample_id               => $sample_id,
                    path                    => $job_id_chain,
                    log                     => $log,
                    sbatch_file_name        => $file_path
                }
            );
        }
    }
    else {

      # Track the number of created xargs scripts per module for Block algorithm
        return
          $xargs_file_counter
          ;
    }
}

1;
