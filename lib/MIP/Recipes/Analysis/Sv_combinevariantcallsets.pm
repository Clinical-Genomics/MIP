package MIP::Recipes::Analysis::Sv_combinevariantcallsets;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_sv_combinevariantcallsets };

}

## Constants
Readonly my $ASTERIX    => q{*};
Readonly my $COLON      => q{:};
Readonly my $DASH       => q{-};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_sv_combinevariantcallsets {

## Function : CombineVariants to combine all structural variants call from different callers.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory {REF}

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
    my $call_type;
    my $family_id;
    my $outaligner_dir;
    my $reference_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        call_type =>
          { default => q{SV}, store => \$call_type, strict_type => 1, },
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
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            store       => \$reference_dir,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_parameters };
    use MIP::Gnu::Coreutils qw(gnu_mv);
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Utility::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_merge bcftools_view bcftools_annotate bcftools_view_and_index_vcf };
    use MIP::Program::Variantcalling::Genmod
      qw{ genmod_annotate genmod_filter };
    use MIP::Program::Variantcalling::Picardtools qw{ sort_vcf };
    use MIP::Program::Variantcalling::Svdb qw{ svdb_merge svdb_query };
    use MIP::Program::Variantcalling::Vcfanno qw{ vcfanno };
    use MIP::Program::Variantcalling::Vt qw{ vt_decompose };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Set::File qw{ set_file_suffix };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Stores the parallel chains that jobIds should be inherited from
    my @parallel_chains;

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            call_type                       => $call_type,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            process_time                    => $time,
            program_directory               => catfile($outaligner_dir),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );
    ## Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      splitpath( $program_info_path . $DOT . q{stderr.txt} );

    ## Assign directories
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir );
    ## Used downstream
    $parameter_href->{$mip_program_name}{indirectory} = $outfamily_directory;

    ## Assign file_tags
    my $outfile_tag =
      $file_info_href->{$family_id}{$mip_program_name}{file_tag};

    ## Paths
    my %file_path_prefix;
    my $merged_file_path_prefix =
      catfile( $temp_directory, $family_id . $UNDERSCORE . $call_type );
    my $outfile_path_prefix =
      catfile( $temp_directory, $family_id . $outfile_tag . $call_type );

    ### Assign suffix
    my %suffix;

    ## Set file suffix for next module within jobid chain
    my $outfile_suffix = set_file_suffix(
        {
            file_suffix => $parameter_href->{$mip_program_name}{outfile_suffix},
            job_id_chain   => $job_id_chain,
            parameter_href => $parameter_href,
            suffix_key     => q{variant_file_suffix},
        }
    );

    ## Collect infiles for all sample_ids for programs that do not do joint calling to enable migration to temporary directory
    _migrate_and_preprocess_single_callers_file(
        {
            active_parameter_href          => $active_parameter_href,
            FILEHANDLE                     => $FILEHANDLE,
            file_info_href                 => $file_info_href,
            file_path_prefix_href          => \%file_path_prefix,
            parallel_chains_ref            => \@parallel_chains,
            parameter_href                 => $parameter_href,
            structural_variant_callers_ref => \@{
                $parameter_href->{dynamic_parameter}{structural_variant_callers}
            },
            suffix_href => \%suffix,
        }
    );

    ## Merged sample files to one family file (samples > 1) else reformat to standardise
    _merge_or_reformat_single_callers_file(
        {
            active_parameter_href          => $active_parameter_href,
            FILEHANDLE                     => $FILEHANDLE,
            file_path_prefix_href          => \%file_path_prefix,
            outfile_suffix                 => $outfile_suffix,
            parameter_href                 => $parameter_href,
            program_info_path              => $program_info_path,
            structural_variant_callers_ref => \@{
                $parameter_href->{dynamic_parameter}{structural_variant_callers}
            },
            suffix_href => \%suffix,
        }
    );

    ## Migrate joint calling per family callers like Manta and Delly
    _migrate_joint_callers_file(
        {
            active_parameter_href          => $active_parameter_href,
            FILEHANDLE                     => $FILEHANDLE,
            file_info_href                 => $file_info_href,
            file_path_prefix_href          => \%file_path_prefix,
            parallel_chains_ref            => \@parallel_chains,
            parameter_href                 => $parameter_href,
            structural_variant_callers_ref => \@{
                $parameter_href->{dynamic_parameter}{structural_variant_callers}
            },
            suffix_href => \%suffix,
        }
    );

    ## Merge structural variant caller's family vcf files
    say {$FILEHANDLE} q{## Merge structural variant caller's family vcf files};

    ## Get parameters
    my @infile_paths;
  STRUCTURAL_CALLER:
    foreach my $structural_variant_caller (
        @{ $parameter_href->{dynamic_parameter}{structural_variant_callers} } )
    {
        ## Expect vcf
        if ( $active_parameter_href->{$structural_variant_caller} ) {

            my $variant_caller_alias =
              $parameter_href->{$structural_variant_caller}{outdir_name};
            push @infile_paths,
              catfile( $temp_directory,
                    $family_id
                  . $UNDERSCORE
                  . $structural_variant_caller
                  . $outfile_suffix
                  . $COLON
                  . $variant_caller_alias );
        }
    }

    svdb_merge(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@infile_paths,
            stdoutfile_path  => $merged_file_path_prefix . $outfile_suffix,
            priority => $active_parameter_href->{sv_svdb_merge_prioritize},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Alternative file tag
    my $alt_file_tag = $EMPTY_STR;

    if ( $active_parameter_href->{sv_vt_decompose} ) {

        ## Update file tag
        $alt_file_tag = $UNDERSCORE . q{vt};

        ## Split multiallelic variants
        say {$FILEHANDLE} q{## Split multiallelic variants};
        vt_decompose(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $merged_file_path_prefix . $outfile_suffix,
                outfile_path => $merged_file_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                smart_decomposition => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    if ( $active_parameter_href->{sv_svdb_query} ) {

        my $infile_path =
          $merged_file_path_prefix . $alt_file_tag . $outfile_suffix;

        ## Update alternative ending
        $alt_file_tag .= $UNDERSCORE . q{svdbq};

        ## Ensure correct infile
        my $annotation_file_counter = 0;

        ## Ensure correct outfiles
        my $outfile_tracker = 0;

        while ( my ( $query_db_file, $query_db_tag ) =
            each %{ $active_parameter_href->{sv_svdb_query_db_files} } )
        {

            if ($annotation_file_counter) {

                $infile_path =
                    $merged_file_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix
                  . $DOT
                  . $outfile_tracker;

                ## Increment now that infile has been set
                $outfile_tracker++;
            }
            svdb_query(
                {
                    bnd_distance    => 25_000,
                    dbfile_path     => $query_db_file,
                    FILEHANDLE      => $FILEHANDLE,
                    frequency_tag   => $query_db_tag . q{AF},
                    hit_tag         => $query_db_tag,
                    infile_path     => $infile_path,
                    stdoutfile_path => $merged_file_path_prefix
                      . $alt_file_tag
                      . $outfile_suffix
                      . $DOT
                      . $outfile_tracker,
                    overlap => 0.8,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
            $annotation_file_counter++;
        }

        ## Rename to remove outfile_tracker
        gnu_mv(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $merged_file_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix
                  . $DOT
                  . $outfile_tracker,
                outfile_path => $merged_file_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Alternative file tag
    my $outfile_alt_file_tag = $alt_file_tag . $UNDERSCORE . q{sorted};

    ## Writes sbatch code to supplied filehandle to sort variants in vcf format
    sort_vcf(
        {
            active_parameter_href => $active_parameter_href,
            FILEHANDLE            => $FILEHANDLE,
            infile_paths_ref =>
              [ $merged_file_path_prefix . $alt_file_tag . $outfile_suffix ],
            outfile => $outfile_path_prefix
              . $outfile_alt_file_tag
              . $outfile_suffix,
            sequence_dict_file => catfile(
                $reference_dir,
                $file_info_href->{human_genome_reference_name_prefix}
                  . $DOT . q{dict}
            ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    $alt_file_tag = $outfile_alt_file_tag;

    ## Remove FILTER ne PASS
    if ( $active_parameter_href->{sv_bcftools_view_filter} ) {

        say {$FILEHANDLE} q{## Remove FILTER ne PASS};
        bcftools_view(
            {
                apply_filters_ref => [qw{ PASS }],
                FILEHANDLE        => $FILEHANDLE,
                infile_path       => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $UNDERSCORE . q{filt}
                  . $outfile_suffix,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Update file tag
        $alt_file_tag .= $UNDERSCORE . q{filt};
    }

    ## Remove common variants
    if ( $active_parameter_href->{sv_genmod_filter} ) {

        my @program_source_commands = get_program_parameters(
            {
                active_parameter_href => $active_parameter_href,
                mip_program_name      => q{genmod},
            }
        );

        say {$FILEHANDLE} q{## Remove common variants};
        genmod_annotate(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => catfile( dirname( devnull() ), q{stdout} ),
                program_source_commands_ref => \@program_source_commands,
                temp_directory_path         => $temp_directory,
                thousand_g_file_path =>
                  $active_parameter_href->{sv_genmod_filter_1000g},
                verbosity => q{v},
            }
        );
        print {$FILEHANDLE} $PIPE . $SPACE;

        ## Update file tag
        $alt_file_tag .= $UNDERSCORE . q{genmod_filter};

        genmod_filter(
            {
                deactive_program_source => 1,
                FILEHANDLE              => $FILEHANDLE,
                infile_path             => $DASH,
                outfile_path            => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                source_main_environment_commands_ref =>
                  \@source_environment_cmds,
                threshold =>
                  $active_parameter_href->{sv_genmod_filter_threshold},
                verbosity => q{v},
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }

    ## Annotate 1000G structural variants
    if ( $active_parameter_href->{sv_vcfanno} ) {

        say {$FILEHANDLE} q{## Annotate 1000G structural variants};
        vcfanno(
            {
                ends        => 1,
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                luafile_path => $active_parameter_href->{sv_vcfanno_lua},
                toml_configfile_path =>
                  $active_parameter_href->{sv_vcfanno_config},
            }
        );
        print {$FILEHANDLE} $PIPE . $SPACE;

        ## Remove "[" and "]" from INFO as it breaks vcf format
        print {$FILEHANDLE}
q?perl -nae 'if($_=~/^#/) {print $_} else {$F[7]=~s/\[||\]//g; print join("\t", @F), "\n"}' ?;

        ## Update file tag
        $alt_file_tag .= $UNDERSCORE . q{vcfanno};

        say {$FILEHANDLE} q{>}
          . $SPACE
          . $outfile_path_prefix
          . $alt_file_tag
          . $outfile_suffix, $NEWLINE;

        if ( $mip_program_mode == 1 ) {

            add_program_outfile_to_sample_info(
                {
                    path             => catfile( $directory, $stderr_file ),
                    program_name     => q{sv_combinevariantcallsets},
                    sample_info_href => $sample_info_href,
                }
            );
        }

        say {$FILEHANDLE}
          q{## Add header for 1000G annotation of structural variants};
        bcftools_annotate(
            {
                FILEHANDLE => $FILEHANDLE,
                headerfile_path =>
                  $active_parameter_href->{sv_vcfannotation_header_lines_file},
                infile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $outfile_suffix,
                outfile_path => $outfile_path_prefix
                  . $alt_file_tag
                  . $UNDERSCORE
                  . q{bcftools_annotate}
                  . $outfile_suffix,
                output_type => q{v},
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ##Update file tag
        $alt_file_tag .= $UNDERSCORE . q{bcftools_annotate};
    }

    ## Then we have something to rename
    if ( $alt_file_tag ne $EMPTY_STR ) {

        ## Writes sbatch code to supplied filehandle to sort variants in vcf format
        sort_vcf(
            {
                active_parameter_href => $active_parameter_href,
                FILEHANDLE            => $FILEHANDLE,
                infile_paths_ref =>
                  [ $outfile_path_prefix . $alt_file_tag . $outfile_suffix ],
                outfile            => $outfile_path_prefix . $outfile_suffix,
                sequence_dict_file => catfile(
                    $reference_dir,
                    $file_info_href->{human_genome_reference_name_prefix}
                      . $DOT . q{dict}
                ),
            }
        );

        say {$FILEHANDLE} $NEWLINE;
    }

    if ( $active_parameter_href->{sv_combinevariantcallsets_bcf_file} ) {

        ## Reformat variant calling file and index
        bcftools_view_and_index_vcf(
            {
                FILEHANDLE          => $FILEHANDLE,
                infile_path         => $outfile_path_prefix . $outfile_suffix,
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => q{b},
            }
        );

        ## Copies file from temporary directory.
        say {$FILEHANDLE} q{## Copy file from temporary directory};
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $outfile_path_prefix . $DOT . q{bcf} . $ASTERIX,
                outfile_path => $outfamily_directory,
            }
        );
    }

    ## Copies file from temporary directory.
    say {$FILEHANDLE} q{## Copy file from temporary directory};
    migrate_file(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $outfile_path_prefix . $outfile_suffix,
            outfile_path => $outfamily_directory,
        }
    );
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $mip_program_mode == 1 ) {

        my $qc_svdb_outfile =
          $family_id . $outfile_tag . $call_type . $outfile_suffix;
        add_program_outfile_to_sample_info(
            {
                path => catfile( $outfamily_directory, $qc_svdb_outfile ),
                program_name     => q{svdb},
                sample_info_href => $sample_info_href,
            }
        );

        $sample_info_href->{sv_vcf_file}{ready_vcf}{path} =
          catfile( $outfamily_directory,
            $family_id . $outfile_tag . $call_type . $outfile_suffix );

        if ( $active_parameter_href->{sv_combinevariantcallsets_bcf_file} ) {

            $sample_info_href->{sv_bcf_file}{path} =
              catfile( $outfamily_directory,
                $family_id . $outfile_tag . $call_type . $DOT . q{bcf} );
        }

        slurm_submit_job_sample_id_dependency_add_to_family(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                parallel_chains_ref     => \@parallel_chains,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

sub _add_to_parallel_chain {

## Function :
## Returns  :
## Arguments: $parallel_chains_ref             => Store structural variant caller parallel chain
##          : $structural_variant_caller_chain => Chain of structural variant caller

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parallel_chains_ref;
    my $structural_variant_caller_chain;

    my $tmpl = {
        parallel_chains_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parallel_chains_ref,
            strict_type => 1,
        },
        structural_variant_caller_chain => {
            defined     => 1,
            required    => 1,
            store       => \$structural_variant_caller_chain,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## If element is not part of array
    if (
        not any {
            $_ eq $structural_variant_caller_chain
        }
        @{$parallel_chains_ref}
      )
    {
        push @{$parallel_chains_ref}, $structural_variant_caller_chain;
    }
    return;
}

sub _migrate_joint_callers_file {

## Function : Migrate joint calling per family callers like Manta and Delly
## Returns  :
## Arguments: $active_parameter_href          => Active parameters for this analysis hash {REF}
##          : $call_type                      => Variant call type
##          : $family_id                      => Family id
##          : $FILEHANDLE                     => Filehandle to write to
##          : $file_info_href                 => File info hash {REF
##          : $file_path_prefix_href          => Store file path prefix {REF}
##          : $outaligner_dir                 => Outaligner_dir used in the analysis
##          : $parallel_chains_ref            => Store structural variant caller parallel chain
##          : $parameter_href                 => Parameter hash {REF}
##          : $structural_variant_callers_ref => Structural variant callers that do not use joint calling
##          : $suffix_href                    => Store suffixes {REF}
##          : $temp_directory                 => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $call_type;
    my $FILEHANDLE;
    my $family_id;
    my $file_info_href;
    my $file_path_prefix_href;
    my $parallel_chains_ref;
    my $parameter_href;
    my $structural_variant_callers_ref;
    my $suffix_href;

    ## Default(s)
    my $outaligner_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        call_type =>
          { default => q{SV}, store => \$call_type, strict_type => 1, },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        FILEHANDLE     => { required => 1, store => \$FILEHANDLE, },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_path_prefix_href,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        parallel_chains_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parallel_chains_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        structural_variant_callers_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$structural_variant_callers_ref,
            strict_type => 1,
        },
        suffix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$suffix_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };

    my $joint_caller = q{pmanta | pdelly_reformat | ptiddit};

  STRUCTURAL_CALLER:
    foreach my $structural_variant_caller ( @{$structural_variant_callers_ref} )
    {

        next STRUCTURAL_CALLER
          if ( not $active_parameter_href->{$structural_variant_caller} );

        next STRUCTURAL_CALLER
          if ( $structural_variant_caller !~ / $joint_caller /xsm );

        ## Expect vcf. Special case: manta, delly, tiddit are processed by joint calling and per family

        ## Assign directories
        my $program_outdirectory_name =
          $parameter_href->{$structural_variant_caller}{outdir_name};
        my $infamily_directory = catfile( $active_parameter_href->{outdata_dir},
            $family_id, $outaligner_dir, $program_outdirectory_name );

        ## Assign file_tags
        my $infile_tag =
          $file_info_href->{$family_id}{$structural_variant_caller}{file_tag};
        my $infile_prefix = $family_id . $infile_tag . $UNDERSCORE . $call_type;

        ## Assign suffix
        my $infile_suffix = get_file_suffix(
            {
                parameter_href => $parameter_href,
                suffix_key     => q{outfile_suffix},
                program_name   => $structural_variant_caller,
            }
        );

        _add_to_parallel_chain(
            {
                parallel_chains_ref => $parallel_chains_ref,
                structural_variant_caller_chain =>
                  $parameter_href->{$structural_variant_caller}{chain},
            }
        );

        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
        migrate_file(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => catfile(
                    $infamily_directory,
                    $infile_prefix . $infile_suffix . $ASTERIX
                ),
                outfile_path => $temp_directory
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        if ( $active_parameter_href->{sv_vt_decompose} ) {

            ## Split multiallelic variants
            say {$FILEHANDLE} q{## Split multiallelic variants};
            vt_decompose(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => catfile(
                        $temp_directory, $infile_prefix . $infile_suffix
                    ),
                    outfile_path => catfile(
                        $temp_directory,
                        $family_id
                          . $UNDERSCORE
                          . $structural_variant_caller
                          . $infile_suffix
                    ),
                    smart_decomposition => 1,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }
    return;
}

sub _migrate_and_preprocess_single_callers_file {

## Function : Collect infiles for all sample_ids for programs that do not do joint calling to enable migration to temporary directory. Add chain of structural variant caller to parallel chains
## Returns  :
## Arguments: $active_parameter_href          => Active parameters for this analysis hash {REF}
##          : $FILEHANDLE                     => Filehandle to write to
##          : $file_info_href                 => File info hash {REF
##          : $file_path_prefix_href          => Store file path prefix {REF}
##          : $outaligner_dir                 => Outaligner_dir used in the analysis
##          : $parallel_chains_ref            => Store structural variant caller parallel chain
##          : $parameter_href                 => Parameter hash {REF}
##          : $structural_variant_callers_ref => Structural variant callers that do not use joint calling
##          : $suffix_href                    => Store suffixes {REF}
##          : $temp_directory                 => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_info_href;
    my $file_path_prefix_href;
    my $parallel_chains_ref;
    my $parameter_href;
    my $structural_variant_callers_ref;
    my $suffix_href;

    ## Default(s)
    my $outaligner_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        FILEHANDLE     => { required => 1, store => \$FILEHANDLE, },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_path_prefix_href,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        parallel_chains_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parallel_chains_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        structural_variant_callers_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$structural_variant_callers_ref,
            strict_type => 1,
        },
        suffix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$suffix_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };

    my $joint_caller = q{pmanta | pdelly_reformat | ptiddit};

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Add merged infile name prefix after merging all BAM files per sample_id
        my $merged_infile_prefix = get_merged_infile_prefix(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

      STRUCTURAL_CALLER:
        foreach
          my $structural_variant_caller ( @{$structural_variant_callers_ref} )
        {

            next STRUCTURAL_CALLER
              if ( not $active_parameter_href->{$structural_variant_caller} );

            next STRUCTURAL_CALLER
              if ( $structural_variant_caller =~ / $joint_caller /xsm );

            ## Expect vcf. Special case: manta, delly and tiddit are processed by joint calling and per family

            ## Assign directories
            my $program_outdirectory_name =
              $parameter_href->{$structural_variant_caller}{outdir_name};
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $outaligner_dir, $program_outdirectory_name );

            ## Assign file_tags
            my $infile_tag =
              $file_info_href->{$sample_id}{$structural_variant_caller}
              {file_tag};
            my $infile_prefix = $merged_infile_prefix . $infile_tag;
            $file_path_prefix_href->{$sample_id}{$structural_variant_caller} =
              catfile( $temp_directory, $infile_prefix );

            ## Assign suffix
            $suffix_href->{$structural_variant_caller} = get_file_suffix(
                {
                    parameter_href => $parameter_href,
                    suffix_key     => q{outfile_suffix},
                    program_name   => $structural_variant_caller,
                }
            );

            _add_to_parallel_chain(
                {
                    parallel_chains_ref => $parallel_chains_ref,
                    structural_variant_caller_chain =>
                      $parameter_href->{$structural_variant_caller}{chain},
                }
            );

            ## Copy file(s) to temporary directory
            say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
            migrate_file(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    infile_path => catfile(
                        $insample_directory,
                        $infile_prefix
                          . $suffix_href->{$structural_variant_caller}
                          . $ASTERIX
                    ),
                    outfile_path => $temp_directory
                }
            );

            say {$FILEHANDLE} q{wait}, $NEWLINE;

            ## Reformat variant calling file and index
            bcftools_view_and_index_vcf(
                {
                    infile_path => $file_path_prefix_href->{$sample_id}
                      {$structural_variant_caller}
                      . $suffix_href->{$structural_variant_caller},
                    outfile_path_prefix =>
                      $file_path_prefix_href->{$sample_id}
                      {$structural_variant_caller},
                    output_type => q{z},
                    FILEHANDLE  => $FILEHANDLE,
                }
            );
        }
    }
    return;
}

sub _merge_or_reformat_single_callers_file {

## Function : Merged sample files to one family file (samples > 1) else reformat to standardise
## Returns  :
## Arguments: $active_parameter_href          => Active parameters for this analysis hash {REF}
##          : $family_id                      => Family ID
##          : $FILEHANDLE                     => Filehandle to write to
##          : $file_path_prefix_href          => Store file path prefix {REF}
##          : $outfile_suffix                 => Outfile suffix
##          : $parameter_href                 => Parameter hash {REF}
##          : $program_info_path              => Program info path
##          : $structural_variant_callers_ref => Structural variant callers that do not use joint calling
##          : $suffix_href                    => Store suffixes {REF}
##          : $temp_directory                 => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $file_path_prefix_href;
    my $outfile_suffix;
    my $parameter_href;
    my $program_info_path;
    my $structural_variant_callers_ref;
    my $suffix_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;

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
        FILEHANDLE            => { required => 1, store => \$FILEHANDLE, },
        file_path_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_path_prefix_href,
            strict_type => 1,
        },
        outfile_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_suffix,
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
            defined     => 1,
            required    => 1,
            store       => \$program_info_path,
            strict_type => 1,
        },
        structural_variant_callers_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$structural_variant_callers_ref,
            strict_type => 1,
        },
        suffix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$suffix_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $joint_caller = q{pmanta | pdelly_reformat | ptiddit};

  STRUCTURAL_CALLER:
    foreach my $structural_variant_caller ( @{$structural_variant_callers_ref} )
    {

        next STRUCTURAL_CALLER
          if ( not $active_parameter_href->{$structural_variant_caller} );

        next STRUCTURAL_CALLER
          if ( $structural_variant_caller =~ / $joint_caller /xsm );

        ## Expect vcf. Special case: joint calling and per family

        ## Assemble file paths by adding file ending
        my @file_paths = map {
                $file_path_prefix_href->{$_}{$structural_variant_caller}
              . $suffix_href->{$structural_variant_caller}
              . $DOT . q{gz}
        } @{ $active_parameter_href->{sample_ids} };

        if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {

            ## Merge all structural variant caller's vcf files per sample_id
            say {$FILEHANDLE}
q{## Merge all structural variant caller's vcf files per sample_id};

            bcftools_merge(
                {
                    FILEHANDLE       => $FILEHANDLE,
                    infile_paths_ref => \@file_paths,
                    outfile_path     => catfile(
                        $temp_directory,
                        $family_id
                          . $UNDERSCORE
                          . $structural_variant_caller
                          . $outfile_suffix
                    ),
                    output_type     => q{v},
                    stderrfile_path => $program_info_path
                      . $UNDERSCORE
                      . $structural_variant_caller
                      . $UNDERSCORE
                      . q{merge.stderr.txt},
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
        else {

            ## Reformat all structural variant caller's vcf files per sample_id
            say {$FILEHANDLE}
q{## Reformat all structural variant caller's vcf files per sample_id};

            bcftools_view(
                {
                    FILEHANDLE   => $FILEHANDLE,
                    infile_path  => $file_paths[0],    #Can be only one
                    outfile_path => catfile(
                        $temp_directory,
                        $family_id
                          . $UNDERSCORE
                          . $structural_variant_caller
                          . $outfile_suffix
                    ),
                    output_type     => q{v},
                    stderrfile_path => $program_info_path
                      . $UNDERSCORE
                      . $structural_variant_caller
                      . $UNDERSCORE
                      . q{merge.stderr.txt},
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }
    return;
}

1;
