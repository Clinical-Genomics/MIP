package MIP::Recipes::Analysis::Cadd;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $CLOSE_PARENTHESIS $COMMA $DASH $DOT $EQUALS $LOG_NAME $NEWLINE $OPEN_PARENTHESIS $PIPE $SEMICOLON $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_cadd analysis_cadd_panel };

}

## Constants
Readonly my $REGION_START  => 2;
Readonly my $REGION_END    => 2;
Readonly my $SEQUENCE_NAME => 1;

sub analysis_cadd {

## Function : Annotate variants with CADD score
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $profile_base_command;
    my $recipe_name;
    my $sample_info_href;

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Bash qw{ gnu_cd gnu_export gnu_unset };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mkdir };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Bcftools qw{ bcftools_annotate bcftools_concat bcftools_view };
    use MIP::Program::Cadd qw{ cadd };
    use MIP::Program::Htslib qw{ htslib_tabix };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };

    my $human_genome_reference_version =
      $file_info_href->{human_genome_reference_version};
    my $assembly_version = _get_cadd_reference_param(
        {
            reference_source  => $file_info_href->{human_genome_reference_source},
            reference_version => $human_genome_reference_version,
        }
    );

    my $cadd_columns_name = join $COMMA, @{ $active_parameter_href->{cadd_column_names} };
    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };
    my $job_id_chain         = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@contigs_size_ordered,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
            }
        )
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };
    my @outfile_paths       = @{ $io{out}{file_paths} };
    my %temp_outdir_path =
      map { $_ => catdir( $active_parameter_href->{temp_directory}, $_ ) }
      @contigs_size_ordered;

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle      = IO::Handle->new();
    my $xargsfilehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### SHELL:
    say {$filehandle} q{## } . $recipe_name;

    ## Add reference dir for CADD mounting point
    my $bash_variable =
        q{MIP_BIND}
      . $EQUALS
      . catdir( $active_parameter_href->{reference_dir},
        qw{ CADD-scripts data annotations} );
    gnu_export(
        {
            bash_variable => $bash_variable,
            filehandle    => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Create temp directories};
  TEMP_OUTDIR:
    foreach my $temp_outdir ( values %temp_outdir_path ) {

        gnu_mkdir(
            {
                filehandle       => $filehandle,
                indirectory_path => $temp_outdir,
            }
        );
        print {$filehandle} $NEWLINE;
    }
    print {$filehandle} $NEWLINE;

    ## View indels and calculate CADD
    say {$filehandle} q{## CADD};

    ## Create file commands for xargs
    my ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number      => $recipe_resource{core_number},
            filehandle       => $filehandle,
            file_path        => $recipe_file_path,
            recipe_info_path => $recipe_info_path,
            xargsfilehandle  => $xargsfilehandle,
        }
    );

    my %cadd_outfile_path;
    ## Process per contig
  CONTIG:
    while ( my ( $index, $contig ) = each @contigs_size_ordered ) {

        ## Get parameters
        $cadd_outfile_path{$contig} = catfile( $temp_outdir_path{$contig},
            $outfile_name_prefix . $DOT . $contig . $DOT . q{tsv.gz} );
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        my $view_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{view} . $DOT . $contig . $outfile_suffix;
        my $view_infile_path = $infile_path{$contig};

        if ( $contig =~ / MT | chrM /xms ) {

            $view_infile_path = $outfile_path_prefix . $DOT . $contig . q{plus.vcf};
            bcftools_concat(
                {
                    filehandle       => $xargsfilehandle,
                    infile_paths_ref => [
                        $infile_path{ $contigs_size_ordered[ $index - 1 ] },
                        $infile_path{$contig}
                    ],
                    outfile_path    => $view_infile_path,
                    output_type     => q{v},
                    rm_dups         => 0,
                    stderrfile_path => $stderrfile_path,
                }
            );
            print {$xargsfilehandle} $SEMICOLON . $SPACE;
        }

        $view_infile_path = _parse_cadd_infile(
            {
                escape_oneliner   => 1,
                filehandle        => $xargsfilehandle,
                infile_path       => $view_infile_path,
                reference_version => $file_info_href->{human_genome_reference_version},
            }
        );

        bcftools_view(
            {
                filehandle      => $xargsfilehandle,
                infile_path     => $view_infile_path,
                types           => q{indels},
                outfile_path    => $view_outfile_path,
                output_type     => q{v},
                stderrfile_path => $stderrfile_path,
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        ## Attempt to escape snakemake lock
        gnu_cd(
            {
                filehandle             => $xargsfilehandle,
                directory_path         => $temp_outdir_path{$contig},
                stderrfile_path_append => $stderrfile_path,
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        cadd(
            {
                filehandle             => $xargsfilehandle,
                genome_build           => $assembly_version,
                infile_path            => $view_outfile_path,
                outfile_path           => $cadd_outfile_path{$contig},
                stderrfile_path_append => $stderrfile_path,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    ### Annotate
    ## Tabix cadd outfile and annotate original vcf file with indel CADD score
    say {$filehandle} q{## Tabix and bcftools annotate};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $recipe_resource{core_number},
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Process per contig
  CONTIG:
    foreach my $contig (@contigs_size_ordered) {

        ## Get parameters
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

        ## Parse outfile in case of grch38
        my $cadd_outfile_path = _parse_cadd_outfile(
            {
                escape_oneliner   => 1,
                filehandle        => $xargsfilehandle,
                infile_path       => $cadd_outfile_path{$contig},
                reference_version => $file_info_href->{human_genome_reference_version},
            }
        );

        ## Create tabix index
        htslib_tabix(
            {
                begin           => $REGION_START,
                end             => $REGION_END,
                filehandle      => $xargsfilehandle,
                force           => 1,
                infile_path     => $cadd_outfile_path,
                sequence        => $SEQUENCE_NAME,
                stderrfile_path => $stderrfile_path,
            }
        );
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        bcftools_annotate(
            {
                annotations_file_path  => $cadd_outfile_path,
                columns_name           => $cadd_columns_name,
                filehandle             => $xargsfilehandle,
                headerfile_path        => $active_parameter_href->{cadd_vcf_header_file},
                infile_path            => $infile_path{$contig},
                outfile_path           => $outfile_path{$contig},
                output_type            => q{v},
                regions_ref            => [$contig],
                stderrfile_path_append => $stderrfile_path,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    gnu_unset(
        {
            bash_variable => q{MIP_BIND},
            filehandle    => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});
    close $xargsfilehandle
      or $log->logcroak(q{Could not close xargsfilehandle});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_paths[0],
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                dependency_method    => q{sample_to_case},
                case_id              => $case_id,
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub analysis_cadd_panel {

## Function : Annotate variants with CADD score for genome build 38
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $profile_base_command;
    my $recipe_name;
    my $sample_info_href;

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Bash qw{ gnu_export gnu_unset };
    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Bcftools qw{ bcftools_annotate bcftools_view };
    use MIP::Program::Cadd qw{ cadd };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

    my $human_genome_reference_version =
      $file_info_href->{human_genome_reference_version};
    my $assembly_version = _get_cadd_reference_param(
        {
            reference_source  => $file_info_href->{human_genome_reference_source},
            reference_version => $human_genome_reference_version,
        }
    );

    my $cadd_columns_name = join $COMMA, @{ $active_parameter_href->{cadd_column_names} };
    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $case_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_path        = $io{out}{file_path};
    my $outfile_suffix      = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### SHELL:
    ## View indels and calculate CADD
    say {$filehandle} q{## CADD};

    ## Add reference dir for CADD mounting point
    my $bash_variable =
        q{MIP_BIND}
      . $EQUALS
      . catdir( $active_parameter_href->{reference_dir},
        qw{ CADD-scripts data annotations} );
    gnu_export(
        {
            bash_variable => $bash_variable,
            filehandle    => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Get parameters
    my $view_infile_path = _parse_cadd_infile(
        {
            filehandle        => $filehandle,
            infile_path       => $infile_path,
            reference_version => $file_info_href->{human_genome_reference_version},
        }
    );
    my $view_outfile_path =
      $outfile_path_prefix . $UNDERSCORE . q{view} . $outfile_suffix;
    bcftools_view(
        {
            filehandle   => $filehandle,
            infile_path  => $view_infile_path,
            types        => q{indels},
            outfile_path => $view_outfile_path,
            output_type  => q{v},
        }
    );
    say {$filehandle} $NEWLINE;

    my $cadd_outfile_path = $outfile_path_prefix . $DOT . q{tsv.gz};
    cadd(
        {
            filehandle   => $filehandle,
            genome_build => $assembly_version,
            infile_path  => $view_outfile_path,
            outfile_path => $cadd_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Parse outfile in case of grch38
    $cadd_outfile_path = _parse_cadd_outfile(
        {
            filehandle        => $filehandle,
            infile_path       => $cadd_outfile_path,
            reference_version => $file_info_href->{human_genome_reference_version},
        }
    );
    say {$filehandle} $NEWLINE;

    ### Annotate
    ## Tabix cadd outfile and annotate original vcf file with indel CADD score
    say {$filehandle} q{## Tabix and bcftools annotate};

    ## Create tabix index
    htslib_tabix(
        {
            begin       => $REGION_START,
            end         => $REGION_END,
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $cadd_outfile_path,
            sequence    => $SEQUENCE_NAME,
        }
    );
    say {$filehandle} $NEWLINE;

    bcftools_annotate(
        {
            annotations_file_path => $cadd_outfile_path,
            columns_name          => $cadd_columns_name,
            filehandle            => $filehandle,
            headerfile_path       => $active_parameter_href->{cadd_vcf_header_file},
            infile_path           => $infile_path,
            outfile_path          => $outfile_path,
            output_type           => q{v},
        }
    );
    say {$filehandle} $NEWLINE;

    gnu_unset(
        {
            bash_variable => q{MIP_BIND},
            filehandle    => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandles
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                dependency_method    => q{sample_to_case},
                case_id              => $case_id,
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _parse_cadd_outfile {

## Function : Parse and return outfile from CADD depending on reference version
## Returns  : $infile_path | $synonyms_file_path
## Arguments: $escape_oneliner   => Escape perl oneliner
##          : $filehandle        => Filehandle
##          : $infile_path       => Cadd outfile path
##          : $reference_version => Genome reference

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $escape_oneliner;
    my $filehandle;
    my $infile_path;
    my $reference_version;

    my $tmpl = {
        escape_oneliner => {
            allow       => [ undef, 0, 1 ],
            store       => \$escape_oneliner,
            strict_type => 1,
        },
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        reference_version => {
            allow       => qr{ \A\d+\z }sxm,
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::File::Path qw{ remove_file_path_suffix };
    use MIP::Program::Htslib qw{ htslib_bgzip };

    return $infile_path if $reference_version eq q{37};

    my $synonyms_file_path = remove_file_path_suffix(
        {
            file_path         => $infile_path,
            file_suffixes_ref => [qw{ .tsv .tsv.gz }],
        }
    );
    $synonyms_file_path .= $UNDERSCORE . q{synonyms.tsv.gz};

    htslib_bgzip(
        {
            decompress      => 1,
            filehandle      => $filehandle,
            force           => 1,
            infile_path     => $infile_path,
            write_to_stdout => 1,
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    perl_nae_oneliners(
        {
            escape_oneliner => $escape_oneliner,
            filehandle      => $filehandle,
            oneliner_name   => q{synonyms_grch37_to_grch38},
        }
    );

    print {$filehandle} $PIPE . $SPACE;

    htslib_bgzip(
        {
            filehandle      => $filehandle,
            force           => 1,
            stdoutfile_path => $synonyms_file_path,
            write_to_stdout => 1,
        }
    );
    print {$filehandle} $SEMICOLON . $SPACE;

    return $synonyms_file_path;
}

sub _parse_cadd_infile {

## Function : Parse CADD infile and preproccess if necessary, depending on reference version
## Returns  : $infile_path | undef
## Arguments: $escape_oneliner   => Escape perl oneliner
##          : $filehandle        => Filehandle
##          : $infile_path       => Infile path
##          : $reference_version => Genome reference

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $escape_oneliner;
    my $filehandle;
    my $infile_path;
    my $reference_version;

    my $tmpl = {
        escape_oneliner => {
            allow       => [ undef, 0, 1 ],
            store       => \$escape_oneliner,
            strict_type => 1,
        },
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        reference_version => {
            allow       => qr{ \A\d+\z }sxm,
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Program::Bcftools qw{ bcftools_view };

    return $infile_path if $reference_version eq q{37};

    bcftools_view(
        {
            filehandle  => $filehandle,
            infile_path => $infile_path,
            output_type => q{v},
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    ## Perl
    my @infile_commands = perl_nae_oneliners(
        {
            escape_oneliner => $escape_oneliner,
            filehandle      => $filehandle,
            oneliner_name   => q{synonyms_grch38_to_grch37},
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    return;
}

sub _get_cadd_reference_param {

## Function : Get the genome assembly and version for CADD
## Returns  : $assembly_version
## Arguments: $reference_source  => Reference source
##          : $reference_version => Reference version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_source;
    my $reference_version;

    my $tmpl = {
        reference_source => {
            defined     => 1,
            required    => 1,
            store       => \$reference_source,
            strict_type => 1,
        },
        reference_version => {
            allow       => qr{ \A\d+\z }sxm,
            store       => \$reference_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %cadd_map = (
        grch37 => q{GRCh37},
        grch38 => q{GRCh38},
        hg19   => q{GRCh37},
    );

    my $assembly_version = $cadd_map{ $reference_source . $reference_version };
    return $assembly_version;
}

1;
