package MIP::Recipes::Analysis::Mip_vcfparser;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_mip_vcfparser analysis_mip_vcfparser_rio };

}

## Constants
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
##          : $infamily_directory      => In family directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
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
    my $infamily_directory;
    my $job_id_href;
    my $outfamily_directory;
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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        file_path          => { strict_type => 1, store => \$file_path },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw(get_core_number);
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Mip_vcfparser qw{ mip_vcfparser };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
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
            process_time                    => $time,
            program_directory               => $outaligner_dir,
            program_name                    => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
            temp_directory                  => $temp_directory,
        }
    );

    # Used downstream
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $infile_tag =
      $file_info_href->{$family_id}{pvarianteffectpredictor}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

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
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
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

    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        my $padding;
        if ( $contig =~ /MT|M/ ) {

            # Special case for mitochondrial contig annotation
            $padding = $MITO_PADDING;
        }

        my @select_feature_annotation_columns;
        my $select_file;
        my $select_file_matching_column;
        my $select_outfile;
        if ( $active_parameter_href->{vcfparser_select_file} ) {

            if (
                !check_entry_hash_of_array(
                    {
                        element  => $contig,
                        hash_ref => $file_info_href,
                        key      => q{select_file_contigs},

                    }
                )
              )
            {

                $select_file =
                  catfile( $active_parameter_href->{vcfparser_select_file} )
                  ;    #List of genes to analyse separately
                $select_file_matching_column = $active_parameter_href
                  ->{vcfparser_select_file_matching_column}
                  ;    #Column of HGNC Symbol in SelectFile (-sf)

                if (
                    (
                        $active_parameter_href
                        ->{vcfparser_select_feature_annotation_columns}
                    )
                    && (
                        @{
                            $active_parameter_href
                              ->{vcfparser_select_feature_annotation_columns}
                        }
                    )
                  )
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
                parse_vep => $active_parameter_href->{pvarianteffectpredictor},
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

    if ( $mip_program_mode == 1 ) {

        ## Clear old vcfparser entry if present
        if ( defined $sample_info_href->{$program_name} ) {

            delete( $sample_info_href->{$program_name} );
        }

        my %gene_panels = (
            range_file  => q{vcfparser_range_feature_file},
            select_file => q{vcfparser_select_file},
        );
        while ( my ( $gene_panel_key, $gene_panel_file ) = each %gene_panels ) {

            ## Collect databases(s) from a potentially merged gene panel file and adds them to sample_info
            collect_gene_panels(
                {
                    aggregate_gene_panel_file =>
                      $active_parameter_href->{$gene_panel_file},
                    aggregate_gene_panels_key => $gene_panel_key,
                    family_id_ref             => $family_id,
                    program_name_ref          => \$program_name,
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
                outdirectory     => $outfamily_directory,
                outfile          => $qc_vcfparser_outfile,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );
    }

    close $XARGSFILEHANDLE;

    my $vcfparser_analysis_type = $EMPTY_STR;
    my @vcfparser_contigs_ref   = \@{ $file_info_href->{contigs_size_ordered} };

    for (
        my $vcfparser_outfile_counter = 0 ;
        $vcfparser_outfile_counter <
        $active_parameter_href->{vcfparser_outfile_count} ;
        $vcfparser_outfile_counter++
      )
    {

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
    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {
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
##          : $infamily_directory      => In family directory
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
##          : $parameter_href          => Parameter hash {REF}
##          : $program_info_path       => The program info path
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path;
    my $infile_lane_prefix_href;
    my $infamily_directory;
    my $job_id_href;
    my $outfamily_directory;
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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
        FILEHANDLE     => { store => \$FILEHANDLE, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        file_path          => { strict_type => 1, store => \$file_path },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        infamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infamily_directory,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw(get_core_number);
    use MIP::Get::File qw{ get_file_suffix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Mip_vcfparser qw{ mip_vcfparser };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Get core number depending on user supplied input exists or not and max number of cores
    $core_number = get_core_number(
        {
            module_core_number   => $core_number,
            modifier_core_number => scalar @{ $file_info_href->{contigs} },
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    # Used downstream
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Tags
    my $infile_tag =
      $file_info_href->{$family_id}{pvarianteffectpredictor}{file_tag};
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

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
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
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

    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        my $padding;
        if ( $contig =~ /MT|M/ ) {

            # Special case for mitochondrial contig annotation
            $padding = $MITO_PADDING;
        }

        my @select_feature_annotation_columns;
        my $select_file;
        my $select_file_matching_column;
        my $select_outfile;
        if ( $active_parameter_href->{vcfparser_select_file} ) {

            if (
                !check_entry_hash_of_array(
                    {
                        element  => $contig,
                        hash_ref => $file_info_href,
                        key      => q{select_file_contigs},

                    }
                )
              )
            {

                $select_file =
                  catfile( $active_parameter_href->{vcfparser_select_file} )
                  ;    #List of genes to analyse separately
                $select_file_matching_column = $active_parameter_href
                  ->{vcfparser_select_file_matching_column}
                  ;    #Column of HGNC Symbol in SelectFile (-sf)

                if (
                    (
                        $active_parameter_href
                        ->{vcfparser_select_feature_annotation_columns}
                    )
                    && (
                        @{
                            $active_parameter_href
                              ->{vcfparser_select_feature_annotation_columns}
                        }
                    )
                  )
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
                parse_vep => $active_parameter_href->{pvarianteffectpredictor},
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

    if ( $mip_program_mode == 1 ) {

        ## Clear old vcfparser entry if present
        if ( defined $sample_info_href->{$program_name} ) {

            delete( $sample_info_href->{$program_name} );
        }

        my %gene_panels = (
            range_file  => q{vcfparser_range_feature_file},
            select_file => q{vcfparser_select_file},
        );
        while ( my ( $gene_panel_key, $gene_panel_file ) = each %gene_panels ) {

            ## Collect databases(s) from a potentially merged gene panel file and adds them to sample_info
            collect_gene_panels(
                {
                    aggregate_gene_panel_file =>
                      $active_parameter_href->{$gene_panel_file},
                    aggregate_gene_panels_key => $gene_panel_key,
                    family_id_ref             => $family_id,
                    program_name_ref          => \$program_name,
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
                outdirectory     => $outfamily_directory,
                outfile          => $qc_vcfparser_outfile,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );
    }

    close $XARGSFILEHANDLE;

    # Track the number of created xargs scripts per module for Block algorithm
    return $xargs_file_counter;
}

1;
