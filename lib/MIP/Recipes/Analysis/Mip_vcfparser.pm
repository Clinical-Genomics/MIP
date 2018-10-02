package MIP::Recipes::Analysis::Mip_vcfparser;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use List::MoreUtils qw { any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ analysis_mip_vcfparser analysis_mip_vcfparser_rio analysis_vcfparser_sv_wes analysis_vcfparser_sv_wgs };

}

## Constants
Readonly my $AMPERSAND    => q{&};
Readonly my $ASTERISK     => q{*};
Readonly my $DOT          => q{.};
Readonly my $EMPTY_STR    => q{};
Readonly my $MITO_PADDING => 10;
Readonly my $NEWLINE      => qq{\n};
Readonly my $SPACE        => q{ };
Readonly my $UNDERSCORE   => q{_};

sub analysis_mip_vcfparser {

## Function : Vcfparser performs parsing of varianteffectpredictor annotated variants
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_info_path;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $call_type;
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        call_type =>
          { default => q{BOTH}, store => \$call_type, strict_type => 1, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path               => { store => \$file_path, strict_type => 1, },
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_info_path =>
          { store => \$program_info_path, strict_type => 1, },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Mip_vcfparser qw{ mip_vcfparser };
    use MIP::QC::Record qw{ add_gene_panel add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    Readonly my $VCFPARSER_OUTFILE_COUNT =>
      $active_parameter_href->{vcfparser_outfile_count} - 1;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number   => $core_number,
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            call_type                       => $call_type,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $outaligner_dir,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_directory = $infamily_directory;

    # Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $infile_tag =
      $file_info_href->{$family_id}{varianteffectpredictor}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            jobid_chain    => $job_id_chain,
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
            job_id_chain   => $job_id_chain,
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
        }
    );

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@{ $file_info_href->{contigs_size_ordered} },
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            indirectory        => $infamily_directory,
            infile             => $infile_prefix,
            program_info_path  => $program_info_path,
            temp_directory     => $temp_directory,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## vcfparser
    say {$FILEHANDLE} q{## vcfparser};

    my $xargs_file_path_prefix;

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        my $padding;
        if ( $contig =~ / MT | M /xms ) {

            # Special case for mitochondrial contig annotation
            $padding = $MITO_PADDING;
        }

        my @select_feature_annotation_columns;
        my $select_file;
        my $select_file_matching_column;
        my $select_outfile;
        if ( $active_parameter_href->{vcfparser_select_file} ) {

            if (
                not check_element_exist_hash_of_array(
                    {
                        element  => $contig,
                        hash_ref => $file_info_href,
                        key      => q{select_file_contigs},

                    }
                )
              )
            {
                # List of genes to analyse separately
                $select_file =
                  catfile( $active_parameter_href->{vcfparser_select_file} );

                # Column of HGNC Symbol in SelectFile (-sf)
                $select_file_matching_column = $active_parameter_href
                  ->{vcfparser_select_file_matching_column};

                if (
                    exists $active_parameter_href
                    ->{vcfparser_select_feature_annotation_columns} )
                {
                    @select_feature_annotation_columns =
                      @{ $active_parameter_href
                          ->{vcfparser_select_feature_annotation_columns} };
                }
                $select_outfile =
                    $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $DOT
                  . q{selected}
                  . $infile_suffix;
            }
        }

        mip_vcfparser(
            {
                FILEHANDLE  => $XARGSFILEHANDLE,
                infile_path => $file_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $infile_suffix,
                padding   => $padding,
                parse_vep => $active_parameter_href->{varianteffectpredictor},
                range_feature_annotation_columns_ref => \@{
                    $active_parameter_href
                      ->{vcfparser_range_feature_annotation_columns}
                },
                range_feature_file_path =>
                  $active_parameter_href->{vcfparser_range_feature_file},
                select_feature_annotation_columns_ref =>
                  \@select_feature_annotation_columns,
                select_feature_file_path       => $select_file,
                select_feature_matching_column => $select_file_matching_column,
                select_outfile                 => $select_outfile,
                stderrfile_path                => $xargs_file_path_prefix
                  . $DOT
                  . $contig
                  . $DOT
                  . q{stderr.txt}
                  . $SPACE,
                stdoutfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $infile_suffix,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## QC Data File(s)
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix,
            outfile_path => $outfamily_directory,

        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    if ( $program_mode == 1 ) {

        ## Clear old vcfparser entry if present
        if ( defined $sample_info_href->{$program_name} ) {

            delete( $sample_info_href->{$program_name} );
        }

        my %gene_panels = (
            range_file  => q{vcfparser_range_feature_file},
            select_file => q{vcfparser_select_file},
        );

      GENE_PANEL:
        while ( my ( $gene_panel_key, $gene_panel_file ) = each %gene_panels ) {

            ## Collect databases(s) from a potentially merged gene panel file and adds them to sample_info
            add_gene_panel(
                {
                    aggregate_gene_panel_file =>
                      $active_parameter_href->{$gene_panel_file},
                    aggregate_gene_panels_key => $gene_panel_key,
                    family_id                 => $family_id,
                    program_name              => $program_name,
                    sample_info_href          => $sample_info_href,
                }
            );
        }

        ## Collect QC metadata info for later use
        my $qc_vcfparser_outfile =
            $outfile_prefix
          . $UNDERSCORE
          . $file_info_href->{contigs_size_ordered}[0]
          . $infile_suffix;
        add_program_outfile_to_sample_info(
            {
                path => catfile( $outfamily_directory, $qc_vcfparser_outfile ),
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );
    }
    close $XARGSFILEHANDLE
      or $log->logcroak(q{Could not close $XARGSFILEHANDLE});

    my $vcfparser_analysis_type = $EMPTY_STR;
    my @vcfparser_contigs_ref   = \@{ $file_info_href->{contigs_size_ordered} };

    ## Determined by vcfparser output
  VCFPARSER_OUTFILE:
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            # SelectFile variants
            $vcfparser_analysis_type = $DOT . q{selected};
            @vcfparser_contigs_ref =
              \@{ $file_info_href->{sorted_select_file_contigs} };
        }

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file(s) from temporary directory};
        ($xargs_file_counter) = xargs_migrate_contig_files(
            {
                contigs_ref => @vcfparser_contigs_ref,
                core_number => $core_number,
                FILEHANDLE  => $FILEHANDLE,
                file_ending => $vcfparser_analysis_type
                  . $infile_suffix
                  . $ASTERISK,
                file_path          => $file_path,
                outdirectory       => $outfamily_directory,
                outfile            => $outfile_prefix,
                program_info_path  => $program_info_path,
                temp_directory     => $temp_directory,
                XARGSFILEHANDLE    => $XARGSFILEHANDLE,
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close $FILEHANDLE});

    if ( $program_mode == 1 ) {
        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

sub analysis_mip_vcfparser_rio {

## Function : Vcfparser performs parsing of varianteffectpredictor annotated variants
## Returns  : |xargs_file_counter
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $FILEHANDLE              => Filehandle to write to
##          : $file_info_href          => File_info hash {REF}
##          : $file_path               => File path
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $stderr_path             => Stderr path of the block script
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_info_path;
    my $program_name;
    my $sample_info_href;
    my $stderr_path;

    ## Default(s)
    my $call_type;
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        call_type => {
            default     => q{BOTH},
            store       => \$call_type,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path => {
            store       => \$file_path,
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_info_path => {
            store       => \$program_info_path,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        stderr_path => {
            defined     => 1,
            required    => 1,
            store       => \$stderr_path,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Hash qw{ check_element_exist_hash_of_array };
    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Mip_vcfparser qw{ mip_vcfparser };
    use MIP::QC::Record qw{ add_gene_panel add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            module_core_number   => $core_number,
        }
    );

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    my $outfamily_directory = $infamily_directory;

    # Used downstream
    $parameter_href->{$program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $infile_tag =
      $file_info_href->{$family_id}{varianteffectpredictor}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$program_name}{file_tag};

    ## Files
    my $infile_prefix  = $family_id . $infile_tag . $call_type;
    my $outfile_prefix = $family_id . $outfile_tag . $call_type;

    ## Paths
    my $file_path_prefix    = catfile( $temp_directory, $infile_prefix );
    my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix    => $parameter_href->{$program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## vcfparser
    say {$FILEHANDLE} q{## vcfparser};

    my $xargs_file_path_prefix;

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        my $padding;
        if ( $contig =~ / MT | M /xms ) {

            # Special case for mitochondrial contig annotation
            $padding = $MITO_PADDING;
        }

        my @select_feature_annotation_columns;
        my $select_file;
        my $select_file_matching_column;
        my $select_outfile;
        if ( $active_parameter_href->{vcfparser_select_file} ) {

            if (
                not check_element_exist_hash_of_array(
                    {
                        element  => $contig,
                        hash_ref => $file_info_href,
                        key      => q{select_file_contigs},

                    }
                )
              )
            {

                # List of genes to analyse separately
                $select_file =
                  catfile( $active_parameter_href->{vcfparser_select_file} );

                # Column of HGNC Symbol in SelectFile (-sf)
                $select_file_matching_column = $active_parameter_href
                  ->{vcfparser_select_file_matching_column};

                if (
                    exists $active_parameter_href
                    ->{vcfparser_select_feature_annotation_columns} )
                {
                    @select_feature_annotation_columns =
                      @{ $active_parameter_href
                          ->{vcfparser_select_feature_annotation_columns} };
                }

                $select_outfile =
                    $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $DOT
                  . q{selected}
                  . $infile_suffix;
            }
        }

        mip_vcfparser(
            {
                FILEHANDLE  => $XARGSFILEHANDLE,
                infile_path => $file_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $infile_suffix,
                padding   => $padding,
                parse_vep => $active_parameter_href->{varianteffectpredictor},
                range_feature_annotation_columns_ref => \@{
                    $active_parameter_href
                      ->{vcfparser_range_feature_annotation_columns}
                },
                range_feature_file_path =>
                  $active_parameter_href->{vcfparser_range_feature_file},
                select_feature_annotation_columns_ref =>
                  \@select_feature_annotation_columns,
                select_feature_file_path       => $select_file,
                select_feature_matching_column => $select_file_matching_column,
                select_outfile                 => $select_outfile,
                stderrfile_path                => $xargs_file_path_prefix
                  . $DOT
                  . $contig
                  . $DOT
                  . q{stderr.txt}
                  . $SPACE,
                stdoutfile_path => $outfile_path_prefix
                  . $UNDERSCORE
                  . $contig
                  . $infile_suffix,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## QC Data File(s)
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path_prefix
              . $UNDERSCORE
              . $file_info_href->{contigs_size_ordered}[0]
              . $infile_suffix,
            outfile_path => $outfamily_directory,

        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    if ( $program_mode == 1 ) {

        ## Clear old vcfparser entry if present
        if ( defined $sample_info_href->{$program_name} ) {

            delete( $sample_info_href->{$program_name} );
        }

        my %gene_panels = (
            range_file  => q{vcfparser_range_feature_file},
            select_file => q{vcfparser_select_file},
        );

      GENE_PANEL:
        while ( my ( $gene_panel_key, $gene_panel_file ) = each %gene_panels ) {

            ## Collect databases(s) from a potentially merged gene panel file and adds them to sample_info
            add_gene_panel(
                {
                    aggregate_gene_panel_file =>
                      $active_parameter_href->{$gene_panel_file},
                    aggregate_gene_panels_key => $gene_panel_key,
                    family_id                 => $family_id,
                    program_name              => $program_name,
                    sample_info_href          => $sample_info_href,
                }
            );
        }

        ## Collect QC metadata info for later use
        my $qc_vcfparser_outfile =
            $outfile_prefix
          . $UNDERSCORE
          . $file_info_href->{contigs_size_ordered}[0]
          . $infile_suffix;
        add_program_outfile_to_sample_info(
            {
                path => catfile( $outfamily_directory, $qc_vcfparser_outfile ),
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );
    }
    close $XARGSFILEHANDLE
      or $log->logcroak(q{Could not close $XARGSFILEHANDLE});

    # Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

sub analysis_vcfparser_sv_wes {

## Function : Vcfparser performs parsing of varianteffectpredictor annotated wes SV variants
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $FILEHANDLE              => Sbatch filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw{ get_core_number };
    use MIP::Get::Analysis qw{ get_vcf_parser_analysis_suffix };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_concatenate_variants };
    use MIP::Program::Variantcalling::Mip_vcfparser qw{ mip_vcfparser };
    use MIP::QC::Record
      qw{ add_most_complete_vcf add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script};

    ### PREPROCESSING:

    ## Constants
    Readonly my $MT_PADDING => 10;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my %io = get_io_files(
        {
            id             => $family_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );

    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path =
      $infile_path_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK;
    my $temp_infile_path_prefix = $io{temp}{file_path_prefix};
    my $temp_infile_path        = $temp_infile_path_prefix . $infile_suffix;

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
            modifier_core_number => 1,
            module_core_number =>
              $active_parameter_href->{module_core_number}{$program_name},
        }
    );

    my @vcfparser_analysis_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count =>
              $active_parameter_href->{sv_vcfparser_outfile_count},
        }
    );

    ## Set and get the io files per chain, id and stream
    my @set_outfile_name_prefixes =
      map { $infile_name_prefix . $_ } @vcfparser_analysis_types;
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $family_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@vcfparser_analysis_types,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                program_name     => $program_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my @outfile_suffixes    = @{ $io{out}{file_suffixes} };

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $infile_path,
            outfile_path => $temp_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    ## vcfparser
    say {$FILEHANDLE} q{## vcfparser};

    my @select_feature_annotation_columns;
    my $select_file;
    my $select_file_matching_column;
    my $select_outfile;
    if ( $active_parameter_href->{sv_vcfparser_select_file} ) {

        ## List of genes to analyse separately
        $select_file =
          catfile( $active_parameter_href->{sv_vcfparser_select_file} );

        ## Column of HGNC Symbol in select file ("-sf")
        $select_file_matching_column =
          $active_parameter_href->{sv_vcfparser_select_file_matching_column};

        if (
            exists $active_parameter_href
            ->{sv_vcfparser_select_feature_annotation_columns} )
        {

            @select_feature_annotation_columns =
              @{ $active_parameter_href
                  ->{sv_vcfparser_select_feature_annotation_columns} };
        }

        ## Select outfile
        $select_outfile = $outfile_path_prefix . $outfile_suffixes[1];
    }

    mip_vcfparser(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $temp_infile_path,
            parse_vep   => $active_parameter_href->{sv_varianteffectpredictor},
            per_gene    => $active_parameter_href->{sv_vcfparser_per_gene},
            range_feature_annotation_columns_ref => \@{
                $active_parameter_href
                  ->{sv_vcfparser_range_feature_annotation_columns}
            },
            range_feature_file_path =>
              $active_parameter_href->{sv_vcfparser_range_feature_file},
            select_feature_annotation_columns_ref =>
              \@select_feature_annotation_columns,
            stderrfile_path => $file_path . $DOT . q{stderr.txt},
            stdoutfile_path => $outfile_path_prefix . $outfile_suffixes[0],
            select_feature_file_path       => $select_file,
            select_feature_matching_column => $select_file_matching_column,
            select_outfile                 => $select_outfile,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
                path             => $select_outfile,
            }
        );

        my %gene_panels = (
            range_file  => q{sv_vcfparser_range_feature_file},
            select_file => q{sv_vcfparser_select_file},
        );
      GENE_PANEL:
        while ( my ( $gene_panel_key, $gene_panel_file ) = each %gene_panels ) {

            ## Collect databases(s) from a potentially merged gene panel file and adds them to sample_info
            add_gene_panel(
                {
                    aggregate_gene_panel_file =>
                      $active_parameter_href->{$gene_panel_file},
                    aggregate_gene_panels_key => $gene_panel_key,
                    family_id => $arg_href->{active_parameter_href}{family_id},
                    program_name     => $program_name,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

sub analysis_vcfparser_sv_wgs {

## Function : Vcfparser performs parsing of varianteffectpredictor annotated ws SV variants
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $FILEHANDLE              => Sbatch filehandle to write to
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Hash qw{ check_element_exist_hash_of_array };
    use MIP::Get::Analysis qw{ get_vcf_parser_analysis_suffix };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_concatenate_variants };
    use MIP::Program::Variantcalling::Mip_vcfparser qw{ mip_vcfparser };
    use MIP::QC::Record
      qw{ add_most_complete_vcf add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script};

    ### PREPROCESSING:

    ## Constants
    Readonly my $MT_PADDING => 10;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    my %io = get_io_files(
        {
            id             => $family_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );

    my $indir_path_prefix  = $io{in}{dir_path_prefix};
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %temp_infile_path   = %{ $io{temp}{file_path_href} };

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};
    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    my $xargs_file_path_prefix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    my @vcfparser_analysis_types = get_vcf_parser_analysis_suffix(
        {
            vcfparser_outfile_count =>
              $active_parameter_href->{sv_vcfparser_outfile_count},
        }
    );

    ## Set and get the io files per chain, id and stream
    my @set_outfile_name_prefixes =
      map { $infile_name_prefix . $_ } @vcfparser_analysis_types;
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $family_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@vcfparser_analysis_types,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                program_name     => $program_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix       = $io{out}{dir_path_prefix};
    my $outfile_path_prefix      = $io{out}{file_path_prefix};
    my $outfile_suffix           = $io{out}{file_suffix};
    my @outfile_suffixes         = @{ $io{out}{file_suffixes} };
    my $temp_outfile_path_prefix = $io{temp}{file_path_prefix};

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    my @contigs = @{ $file_info_href->{contigs} };

    ## Copy file(s) to temporary directory
    say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    ($xargs_file_counter) = xargs_migrate_contig_files(
        {
            contigs_ref        => \@contigs,
            core_number        => $core_number,
            indirectory        => $indir_path_prefix,
            infile             => $infile_name_prefix,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
            temp_directory     => $temp_directory,
        }
    );

    ## vcfparser
    say {$FILEHANDLE} q{## vcfparser};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig (@contigs) {

        my $padding;
        if ( $contig =~ / MT | M /sxm ) {

            ## Special case for mitochondrial contig annotation
            $padding = $MT_PADDING;
        }

        my $vcfparser_temp_outfile_path =
          $temp_outfile_path_prefix . $DOT . $contig . $outfile_suffix;
        my $vcfparser_xargs_file_path_prefix =
          $xargs_file_path_prefix . $DOT . $contig;
        my @select_feature_annotation_columns;
        my $select_file;
        my $select_file_matching_column;
        my $select_outfile;
        if ( $active_parameter_href->{sv_vcfparser_select_file} ) {

            if (
                not check_element_exist_hash_of_array(
                    {
                        element  => $contig,
                        hash_ref => $file_info_href,
                        key      => q{select_file_contigs},
                    }
                )
              )
            {

                ## List of genes to analyse separately
                $select_file =
                  catfile( $active_parameter_href->{sv_vcfparser_select_file} );

                ## Column of HGNC Symbol in select file ("-sf")
                $select_file_matching_column = $active_parameter_href
                  ->{sv_vcfparser_select_file_matching_column};

                if (
                    exists $active_parameter_href
                    ->{sv_vcfparser_select_feature_annotation_columns} )
                {

                    @select_feature_annotation_columns =
                      @{ $active_parameter_href
                          ->{sv_vcfparser_select_feature_annotation_columns} };
                }

                ## Select outfile
                $select_outfile =
                    $temp_outfile_path_prefix
                  . $DOT
                  . $contig
                  . $outfile_suffixes[1];
            }
        }

        mip_vcfparser(
            {
                FILEHANDLE  => $XARGSFILEHANDLE,
                infile_path => $temp_infile_path{$contig},
                padding     => $padding,
                parse_vep =>
                  $active_parameter_href->{sv_varianteffectpredictor},
                per_gene => $active_parameter_href->{sv_vcfparser_per_gene},
                range_feature_annotation_columns_ref => \@{
                    $active_parameter_href
                      ->{sv_vcfparser_range_feature_annotation_columns}
                },
                range_feature_file_path =>
                  $active_parameter_href->{sv_vcfparser_range_feature_file},
                select_feature_annotation_columns_ref =>
                  \@select_feature_annotation_columns,
                stderrfile_path => $vcfparser_xargs_file_path_prefix
                  . $DOT
                  . q{stderr.txt},
                stdoutfile_path                => $vcfparser_temp_outfile_path,
                select_feature_file_path       => $select_file,
                select_feature_matching_column => $select_file_matching_column,
                select_outfile                 => $select_outfile,
            }
        );
        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## QC Data File(s)
    migrate_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $temp_outfile_path_prefix
              . $DOT
              . $contigs[0]
              . $outfile_suffix
              . $ASTERISK,
            outfile_path => $outdir_path_prefix,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

  ANALYSIS_SUFFIXES:
    foreach my $analysis_suffix (@outfile_suffixes) {
        gatk_concatenate_variants(
            {
                active_parameter_href => $active_parameter_href,
                elements_ref          => \@contigs,
                continue              => 1,
                FILEHANDLE            => $FILEHANDLE,
                infile_postfix        => $analysis_suffix,
                infile_prefix         => $temp_outfile_path_prefix,
                outfile_path_prefix   => $outfile_path_prefix,
                outfile_suffix        => $analysis_suffix,
            }
        );
    }
    say {$FILEHANDLE} q{wait} . $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    close $XARGSFILEHANDLE
      or $log->logcroak(q{Could not close XARGSFILEHANDLE});

    if ( $program_mode == 1 ) {

        my $outfile_sample_info_prefix =
          $outfile_path_prefix . $DOT . $contigs[0] . $outfile_suffix;

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
                path             => $outfile_sample_info_prefix,
            }
        );

        my %gene_panels = (
            range_file  => q{sv_vcfparser_range_feature_file},
            select_file => q{sv_vcfparser_select_file},
        );
      GENE_PANEL:
        while ( my ( $gene_panel_key, $gene_panel_file ) = each %gene_panels ) {

            ## Collect databases(s) from a potentially merged gene panel file and adds them to sample_info
            add_gene_panel(
                {
                    aggregate_gene_panel_file =>
                      $active_parameter_href->{$gene_panel_file},
                    aggregate_gene_panels_key => $gene_panel_key,
                    family_id => $arg_href->{active_parameter_href}{family_id},
                    program_name     => $program_name,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

1;
