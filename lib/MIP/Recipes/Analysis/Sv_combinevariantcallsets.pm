package MIP::Recipes::Analysis::Sv_combinevariantcallsets;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $COLON $DOT $EMPTY_STR $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_sv_combinevariantcallsets };

}

sub analysis_sv_combinevariantcallsets {

## Function : CombineVariants to combine all structural variants call from different callers
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools
      qw{ bcftools_merge bcftools_norm bcftools_view bcftools_view_and_index_vcf };
    use MIP::Program::Svdb qw{ svdb_merge };
    use MIP::Sample_info qw{ set_file_path_to_store
      set_recipe_outfile_in_sample_info
      set_recipe_metafile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Stores the parallel chains that job ids should be inherited from
    my @parallel_chains;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my @structural_variant_callers;

    ## Only process active callers
    foreach
      my $structural_variant_caller ( @{ $parameter_href->{cache}{structural_variant_callers} } )
    {
        if ( $active_parameter_href->{$structural_variant_caller} ) {

            push @structural_variant_callers, $structural_variant_caller;
        }
    }

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $recipe{job_id_chain},
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
        }
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path        = $io{out}{file_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );
    ## Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      splitpath( $recipe_info_path . $DOT . q{stderr.txt} );

    ### SHELL:

    ## Collect infiles for all sample_ids for programs that do not do joint calling to enable migration to temporary directory
    # Paths for structural variant callers to be merged
    my %file_path;

    _preprocess_single_callers_file(
        {
            active_parameter_href          => $active_parameter_href,
            filehandle                     => $filehandle,
            file_info_href                 => $file_info_href,
            file_path_href                 => \%file_path,
            parallel_chains_ref            => \@parallel_chains,
            parameter_href                 => $parameter_href,
            structural_variant_callers_ref => \@structural_variant_callers,
        }
    );

    ## Merged sample files to one case file (samples > 1) else reformat to standardise
    _merge_or_reformat_single_callers_file(
        {
            active_parameter_href          => $active_parameter_href,
            filehandle                     => $filehandle,
            file_path_href                 => \%file_path,
            outdir_path_prefix             => $outdir_path_prefix,
            outfile_suffix                 => $outfile_suffix,
            parameter_href                 => $parameter_href,
            recipe_info_path               => $recipe_info_path,
            structural_variant_callers_ref => \@structural_variant_callers,
        }
    );

    ## Migrate joint calling per case callers like Manta and Delly
    _preprocess_joint_callers_file(
        {
            active_parameter_href          => $active_parameter_href,
            filehandle                     => $filehandle,
            file_info_href                 => $file_info_href,
            file_path_href                 => \%file_path,
            outdir_path_prefix             => $outdir_path_prefix,
            outfile_suffix                 => $outfile_suffix,
            parallel_chains_ref            => \@parallel_chains,
            parameter_href                 => $parameter_href,
            structural_variant_callers_ref => \@structural_variant_callers,
        }
    );

    ## Merge structural variant caller's case vcf files
    say {$filehandle} q{## Merge structural variant caller's case vcf files};

    ## Get parameters
    my @svdb_infile_paths;
  STRUCTURAL_CALLER:
    foreach my $structural_variant_caller (@structural_variant_callers) {

        ## Only use first part of name
        my ($variant_caller_prio_tag) = split /_/sxm, $structural_variant_caller;
        push @svdb_infile_paths,
          catfile( $outdir_path_prefix,
                $case_id
              . $UNDERSCORE
              . $structural_variant_caller
              . $outfile_suffix
              . $COLON
              . $variant_caller_prio_tag );
    }

    svdb_merge(
        {
            filehandle       => $filehandle,
            infile_paths_ref => \@svdb_infile_paths,
            priority         => $active_parameter_href->{sv_svdb_merge_prioritize},
            same_order       => 1,
            stdoutfile_path  => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Alternative file tag
    my $alt_file_tag = $EMPTY_STR;

    if ( $active_parameter_href->{sv_decompose} ) {

        ## Update file tag
        $alt_file_tag = $UNDERSCORE . q{decompose};

        ## Split multiallelic variants
        say {$filehandle} q{## Split multiallelic variants};
        bcftools_norm(
            {
                filehandle   => $filehandle,
                infile_path  => $outfile_path,
                multiallelic => q{-},
                outfile_path => $outfile_path_prefix . $alt_file_tag . $outfile_suffix,
            }
        );
        say {$filehandle} $NEWLINE;

        gnu_mv(
            {
                filehandle   => $filehandle,
                force        => 1,
                infile_path  => $outfile_path_prefix . $alt_file_tag . $outfile_suffix,
                outfile_path => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;

    }

    if ( $active_parameter_href->{sv_combinevariantcallsets_bcf_file} ) {

        ## Reformat variant calling file and index
        bcftools_view_and_index_vcf(
            {
                filehandle          => $filehandle,
                infile_path         => $outfile_path,
                index_type          => q{csi},
                outfile_path_prefix => $outfile_path_prefix,
                output_type         => q{b},
            }
        );
    }

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        $sample_info_href->{sv_vcf_file}{ready_vcf}{path} = $outfile_path;

        if ( $active_parameter_href->{sv_combinevariantcallsets_bcf_file} ) {

            my $sv_bcf_file_path = $outfile_path_prefix . $DOT . q{bcf};
            set_recipe_metafile_in_sample_info(
                {
                    metafile_tag     => q{sv_bcf_file},
                    path             => $sv_bcf_file_path,
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );

            set_file_path_to_store(
                {
                    format           => q{bcf},
                    id               => $case_id,
                    path             => $sv_bcf_file_path,
                    path_index       => $sv_bcf_file_path . $DOT . q{csi},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                parallel_chains_ref => \@parallel_chains,
                recipe_file_path    => $recipe_file_path,
                sample_ids_ref      => \@{ $active_parameter_href->{sample_ids} },
                submission_profile  => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _add_to_parallel_chain {

## Function : Add to parallel chain
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

sub _preprocess_joint_callers_file {

## Function : Preprocess joint calling per case callers like Manta and Delly. And store merged outfile per caller
## Returns  :
## Arguments: $active_parameter_href          => Active parameters for this analysis hash {REF}
##          : $case_id                        => Family id
##          : $filehandle                     => Filehandle to write to
##          : $file_info_href                 => File info hash {REF
##          : $file_path_href                 => Store file path prefix {REF}
##          : $outdir_path_prefix             => Outdir path prefix
##          : $outfile_suffix                 => Outfile suffix
##          : $parallel_chains_ref            => Store structural variant caller parallel chain
##          : $parameter_href                 => Parameter hash {REF}
##          : $structural_variant_callers_ref => Structural variant callers that do not use joint calling

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $filehandle;
    my $file_info_href;
    my $file_path_href;
    my $outdir_path_prefix;
    my $outfile_suffix;
    my $parallel_chains_ref;
    my $parameter_href;
    my $structural_variant_callers_ref;

    ## Default(s)
    my $case_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        filehandle     => { required => 1, store => \$filehandle, },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_path_href,
            strict_type => 1,
        },
        outdir_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path_prefix,
            strict_type => 1,
        },
        outfile_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_suffix,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files };

    my $joint_caller = q{manta | delly_reformat | tiddit};
    my $stream       = q{out};

  STRUCTURAL_CALLER:
    foreach my $structural_variant_caller ( @{$structural_variant_callers_ref} ) {

        next STRUCTURAL_CALLER
          if ( $structural_variant_caller !~ / $joint_caller /xsm );

        ## Expect vcf. Special case: manta, delly, tiddit are processed by joint calling and per case

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $case_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $structural_variant_caller,
                stream         => $stream,
            }
        );

        my $infile_path_prefix = $sample_io{$stream}{file_path_prefix};
        my $infile_suffix      = $sample_io{$stream}{file_suffix};
        my $infile_path        = $infile_path_prefix . $infile_suffix;

        _add_to_parallel_chain(
            {
                parallel_chains_ref             => $parallel_chains_ref,
                structural_variant_caller_chain =>
                  $parameter_href->{$structural_variant_caller}{chain},
            }
        );

        my $decompose_outfile_path = catfile( $outdir_path_prefix,
            $case_id . $UNDERSCORE . $structural_variant_caller . $outfile_suffix );
        ## Store merged outfile per caller
        push @{ $file_path_href->{$structural_variant_caller} }, $decompose_outfile_path;

        if ( $active_parameter_href->{sv_decompose} ) {

            ## Split multiallelic variants
            say {$filehandle} q{## Split multiallelic variants};
            bcftools_norm(
                {
                    filehandle   => $filehandle,
                    infile_path  => $infile_path,
                    multiallelic => q{-},
                    outfile_path => $decompose_outfile_path,
                }
            );
            say {$filehandle} $NEWLINE;
        }
    }
    return;
}

sub _preprocess_single_callers_file {

## Function : Collect infiles for all sample_ids for programs that do not do joint calling. Add chain of structural variant caller to parallel chains
## Returns  :
## Arguments: $active_parameter_href          => Active parameters for this analysis hash {REF}
##          : $filehandle                     => Filehandle to write to
##          : $file_info_href                 => File info hash {REF
##          : $file_path_href                 => Store file path prefix {REF}
##          : $parallel_chains_ref            => Store structural variant caller parallel chain
##          : $parameter_href                 => Parameter hash {REF}
##          : $structural_variant_callers_ref => Structural variant callers that do not use joint calling
##          : $temp_directory                 => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $filehandle;
    my $file_info_href;
    my $file_path_href;
    my $parallel_chains_ref;
    my $parameter_href;
    my $structural_variant_callers_ref;

    ## Default(s)
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        filehandle     => { required => 1, store => \$filehandle, },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_path_href,
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files };

    my $joint_caller = q{manta | delly_reformat | tiddit};
    my $stream       = q{out};

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

      STRUCTURAL_CALLER:
        foreach my $structural_variant_caller ( @{$structural_variant_callers_ref} ) {

            next STRUCTURAL_CALLER
              if ( $structural_variant_caller =~ / $joint_caller /xsm );

            ## Expect vcf. Special case: manta, delly and tiddit are processed by joint calling and per case

            ## Get the io infiles per chain and id
            my %sample_io = get_io_files(
                {
                    id             => $sample_id,
                    file_info_href => $file_info_href,
                    parameter_href => $parameter_href,
                    recipe_name    => $structural_variant_caller,
                    stream         => $stream,
                }
            );

            my $infile_path_prefix = $sample_io{$stream}{file_path_prefix};
            my $infile_suffix      = $sample_io{$stream}{file_suffix};
            my $infile_path        = $infile_path_prefix . $infile_suffix;

            push @{ $file_path_href->{$structural_variant_caller} }, $infile_path . $DOT . q{gz};

            _add_to_parallel_chain(
                {
                    parallel_chains_ref             => $parallel_chains_ref,
                    structural_variant_caller_chain =>
                      $parameter_href->{$structural_variant_caller}{chain},
                }
            );

            ## Reformat variant calling file and index
            bcftools_view_and_index_vcf(
                {
                    infile_path         => $infile_path,
                    outfile_path_prefix => $infile_path_prefix,
                    output_type         => q{z},
                    filehandle          => $filehandle,
                }
            );
        }
    }
    return;
}

sub _merge_or_reformat_single_callers_file {

## Function : Merged sample files to one case file (samples > 1) else reformat to standardise
## Returns  :
## Arguments: $active_parameter_href          => Active parameters for this analysis hash {REF}
##          : $case_id                        => Family ID
##          : $filehandle                     => Filehandle to write to
##          : $file_path_href                 => Store file path prefix {REF}
##          : $outdir_path_prefix             => Outdir path prefix
##          : $outfile_suffix                 => Outfile suffix
##          : $parameter_href                 => Parameter hash {REF}
##          : $recipe_info_path               => Program info path
##          : $structural_variant_callers_ref => Structural variant callers that do not use joint calling

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $filehandle;
    my $file_path_href;
    my $outdir_path_prefix;
    my $outfile_suffix;
    my $parameter_href;
    my $recipe_info_path;
    my $structural_variant_callers_ref;

    ## Default(s)
    my $case_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        filehandle     => { required => 1, store => \$filehandle, },
        file_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_path_href,
            strict_type => 1,
        },
        outdir_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path_prefix,
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
        recipe_info_path => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_info_path,
            strict_type => 1,
        },
        structural_variant_callers_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$structural_variant_callers_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $joint_caller = q{manta | delly_reformat | tiddit};

  STRUCTURAL_CALLER:
    foreach my $structural_variant_caller ( @{$structural_variant_callers_ref} ) {

        next STRUCTURAL_CALLER
          if ( $structural_variant_caller =~ / $joint_caller /xsm );

        ## Expect vcf. Special case: joint calling and per case

        ## Assemble file paths by adding file ending
        my @merge_infile_paths = @{ $file_path_href->{$structural_variant_caller} };
        my $merge_outfile_path = catfile( $outdir_path_prefix,
            $case_id . $UNDERSCORE . $structural_variant_caller . $outfile_suffix );
        ## Store merged outfile per caller
        push @{ $file_path_href->{$structural_variant_caller} }, $merge_outfile_path;

        if ( scalar @{ $active_parameter_href->{sample_ids} } > 1 ) {

            ## Merge all structural variant caller's vcf files per sample_id
            say {$filehandle} q{## Merge all structural variant caller's vcf files per sample_id};

            bcftools_merge(
                {
                    filehandle       => $filehandle,
                    infile_paths_ref => \@merge_infile_paths,
                    outfile_path     => $merge_outfile_path,
                    output_type      => q{v},
                    stderrfile_path  => $recipe_info_path
                      . $UNDERSCORE
                      . $structural_variant_caller
                      . $UNDERSCORE
                      . q{merge.stderr.txt},
                }
            );
            say {$filehandle} $NEWLINE;
        }
        else {

            ## Reformat all structural variant caller's vcf files per sample_id
            say {$filehandle}
              q{## Reformat all structural variant caller's vcf files per sample_id};

            bcftools_view(
                {
                    filehandle      => $filehandle,
                    infile_path     => $merge_infile_paths[0],    # Can be only one
                    outfile_path    => $merge_outfile_path,
                    output_type     => q{v},
                    stderrfile_path => $recipe_info_path
                      . $UNDERSCORE
                      . $structural_variant_caller
                      . $UNDERSCORE
                      . q{merge.stderr.txt},
                }
            );
            say {$filehandle} $NEWLINE;
        }
    }
    return;
}

1;
