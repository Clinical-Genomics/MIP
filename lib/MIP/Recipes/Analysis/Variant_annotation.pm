package MIP::Recipes::Analysis::Variant_annotation;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{
  $AMPERSAND
  $ASTERISK
  $CLOSE_BRACKET
  $CLOSE_PARENTHESIS
  $DASH
  $DOLLAR_SIGN
  $DOT
  $EMPTY_STR
  $FORWARD_SLASH
  $LOG_NAME
  $NEWLINE
  $OPEN_BRACKET
  $OPEN_PARENTHESIS
  $PIPE
  $SPACE
  $UNDERSCORE
};

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      analysis_variant_annotation
      analysis_variant_annotation_panel
      analysis_vcfanno_preop
    };

}

sub analysis_variant_annotation {

## Function : Annotate vcf with allelle frequencies
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
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;

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
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Bcftools qw{ bcftools_concat bcftools_index bcftools_view };
    use MIP::Program::Vcfanno qw{ vcfanno };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
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
    my %infile_path        = %{ $io{in}{file_path_href} };

    my @contigs_size_ordered = @{ $file_info_href->{contigs_size_ordered} };
    my $job_id_chain         = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
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

    my @outfile_paths       = @{ $io{out}{file_paths} };
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };
    my $outfile_suffix      = $io{out}{file_suffix};

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

  CONTIG:
    foreach my $contig (@contigs_size_ordered) {

        ## Get parameters
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

        vcfanno(
            {
                filehandle           => $xargsfilehandle,
                infile_path          => $infile_path{$contig},
                luafile_path         => $active_parameter_href->{vcfanno_functions},
                stderrfile_path      => $stderrfile_path,
                toml_configfile_path => $active_parameter_href->{vcfanno_config},
            }
        );

        print {$xargsfilehandle} $PIPE . $SPACE;

        bcftools_view(
            {
                filehandle             => $xargsfilehandle,
                infile_path            => $DASH,
                outfile_path           => $outfile_path{$contig},
                output_type            => q{z},
                stderrfile_path_append => $stderrfile_path,
            }
        );
        say {$xargsfilehandle} $NEWLINE;

    }

    say {$filehandle} q{## Index outfiles};

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

  CONTIG:
    foreach my $contig (@contigs_size_ordered) {
        ## Get parameters
        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

        bcftools_index(
            {
                filehandle      => $xargsfilehandle,
                infile_path     => $outfile_path{$contig},
                output_type     => q{tbi},
                stderrfile_path => $stderrfile_path,
            }
        );
        say {$xargsfilehandle} $NEWLINE;
    }

    close $xargsfilehandle or $log->logcroak(q{Could not close xargsfilehandle});

    say {$filehandle} q{## Concatenate outfiles};

    bcftools_concat(
        {
            filehandle       => $filehandle,
            infile_paths_ref => \@outfile_paths,
            output_type      => q{z},
            outfile_path     => $outfile_path_prefix . $DOT . q{vcf.gz},
            rm_dups          => 0,
            threads          => $recipe_resource{core_number} - 1,
        }
    );
    say {$filehandle} $NEWLINE;

    bcftools_index(
        {
            infile_path => $outfile_path_prefix . $DOT . q{vcf.gz},
            filehandle  => $filehandle,
            output_type => q{tbi},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path_prefix . $DOT . q{vcf.gz},
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{sample_to_case},
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

sub analysis_variant_annotation_panel {

## Function : Annotate vcf with allelle frequencies
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
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;

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
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Program::Bcftools qw{ bcftools_index bcftools_view };
    use MIP::Program::Vcfanno qw{ vcfanno };
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

    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
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
    my $outfile_path = $io{out}{file_path};

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

    say {$filehandle} q{## } . $recipe_name;

    vcfanno(
        {
            filehandle           => $filehandle,
            infile_path          => $infile_path,
            luafile_path         => $active_parameter_href->{vcfanno_functions},
            toml_configfile_path => $active_parameter_href->{vcfanno_config},
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    bcftools_view(
        {
            filehandle   => $filehandle,
            infile_path  => $DASH,
            outfile_path => $outfile_path,
            output_type  => q{z},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Index outfiles};

    bcftools_index(
        {
            filehandle  => $filehandle,
            infile_path => $outfile_path,
            output_type => q{tbi},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandleS
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
                case_id              => $case_id,
                dependency_method    => q{sample_to_case},
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

sub analysis_vcfanno_preop {

## Function: Split annotation file on available cores and launch vcfanno process for each split
## Returns:
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $annotation_href       => Annotation in toml format {REF}
##          : $case_id               => Family id
##          : $infile_path           => Infile
##          : $job_id_href           => Job id hash {REF}
##          : $outfile_path          => Outfile
##          : $parameter_href        => Parameter hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Recipe name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $annotation_href;
    my $infile_path;
    my $job_id_href;
    my $parameter_href;

    ## Default(s)
    my $case_id;
    my $outfile_path;
    my $profile_base_command;
    my $recipe_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        annotation_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$annotation_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        outfile_path => {
            default     => $arg_href->{infile_path},
            store       => \$outfile_path,
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
        recipe_name =>
          { default => q{vcfanno_preops}, store => \$recipe_name, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Io::Write qw{ write_to_file };
    use MIP::Parse::File qw{ parse_file_suffix };
    use MIP::Program::Bcftools qw{ bcftools_concat bcftools_index bcftools_view };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cat gnu_rm gnu_split };
    use MIP::Program::Vcfanno qw{ vcfanno };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

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
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### SHELL:

    ## Set file paths
    my ( $infile_name, $reference_dir ) = fileparse($infile_path);
    my $outdir_path =
      catdir( $active_parameter_href->{outdata_dir}, $case_id, $recipe_name );
    my $temp_outfile_path        = catfile( $outdir_path, $infile_name );
    my $temp_outfile_path_prefix = parse_file_suffix(
        {
            file_name   => $temp_outfile_path,
            file_suffix => q{.vcf},
        }
    );
    my $max_file_number = $recipe_resource{core_number} - 1;

    my @splitted_infiles = _build_file_splitted_file_paths(
        {
            file_path_prefix => $temp_outfile_path_prefix,
            max_file_number  => $max_file_number,
        }
    );

    my @splitted_outfiles = _build_file_splitted_file_paths(
        {
            file_path_prefix => $temp_outfile_path_prefix,
            file_suffix      => q{.vcf.gz},
            max_file_number  => $max_file_number,
        }
    );

    my @splitted_stderrfiles = _build_file_splitted_file_paths(
        {
            file_path_prefix => $recipe_info_path,
            file_suffix      => q{.stderr.txt},
            max_file_number  => $max_file_number,
        }
    );

    say {$filehandle} q{## } . $recipe_name;

    ## Write vcfanno config
    my $vcfanno_config_path = catfile( $outdir_path,
        $recipe_name . $UNDERSCORE . q{config} . $DOT . $random_integer . q{.toml} );
    write_to_file(
        {
            data_href => $annotation_href,
            format    => q{toml},
            path      => $vcfanno_config_path,
        }
    );

    ## Separate header and vcf body
    my $vcf_header_path = $temp_outfile_path_prefix . q{_header.txt};
    my $vcf_body_path   = $temp_outfile_path_prefix . q{_body.txt};
    bcftools_view(
        {
            filehandle   => $filehandle,
            header_only  => 1,
            infile_path  => $infile_path,
            outfile_path => $vcf_header_path,
        }
    );
    print {$filehandle} $NEWLINE;
    bcftools_view(
        {
            filehandle   => $filehandle,
            infile_path  => $infile_path,
            no_header    => 1,
            outfile_path => $vcf_body_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Split Vcf body into the same number as available cores};
    gnu_split(
        {
            filehandle       => $filehandle,
            files            => $recipe_resource{core_number},
            infile_path      => $vcf_body_path,
            numeric_suffixes => 1,
            prefix           => $temp_outfile_path_prefix . $DOT,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Annotating splitted files};
  FILE_NUMBER:
    for my $file_number ( 0 .. $max_file_number ) {

        gnu_cat(
            {
                filehandle       => $filehandle,
                infile_paths_ref => [ $vcf_header_path, $splitted_infiles[$file_number] ],

            }
        );
        print {$filehandle} $PIPE . $SPACE;

        vcfanno(
            {
                filehandle           => $filehandle,
                infile_path          => catfile( $FORWARD_SLASH, qw{ dev stdin } ),
                luafile_path         => $active_parameter_href->{vcfanno_functions},
                toml_configfile_path => $vcfanno_config_path,
                stderrfile_path      => $splitted_stderrfiles[$file_number],
            }
        );
        print {$filehandle} $PIPE . $SPACE;

        bcftools_view(
            {
                filehandle             => $filehandle,
                infile_path            => $DASH,
                outfile_path           => $splitted_outfiles[$file_number],
                output_type            => q{z},
                stderrfile_path_append => $splitted_stderrfiles[$file_number],
            }
        );
        print {$filehandle} $AMPERSAND . $NEWLINE;
    }
    say {$filehandle} q{wait} . $NEWLINE;

    say {$filehandle} q{## Indexing splitted files};
  FILE_NUMBER:
    for my $file_number ( 0 .. $max_file_number ) {

        bcftools_index(
            {
                filehandle             => $filehandle,
                infile_path            => $splitted_outfiles[$file_number],
                output_type            => q{tbi},
                stderrfile_path_append => $splitted_stderrfiles[$file_number],
            }
        );
        print {$filehandle} $AMPERSAND . $NEWLINE;
    }
    say {$filehandle} q{wait} . $NEWLINE;

    say {$filehandle} q{## Concatenate splitted files};
    bcftools_concat(
        {
            filehandle       => $filehandle,
            infile_paths_ref => \@splitted_outfiles,
            outfile_path     => $temp_outfile_path,
            output_type      => q{z},
            rm_dups          => 0,
            threads          => $recipe_resource{core_number},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Index outfile};
    bcftools_index(
        {
            filehandle  => $filehandle,
            infile_path => $temp_outfile_path,
            output_type => q{tbi},
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Check outfile and move};
    _check_outfile_and_move(
        {
            filehandle        => $filehandle,
            outfile_path      => $outfile_path,
            reference_dir     => $reference_dir,
            temp_outfile_path => $temp_outfile_path,
            vcf_id_tag        => $annotation_href->{annotation}[0]{names}[0],
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Cleanup};
    my $splitted_files = join $SPACE,
      ( @splitted_infiles, @splitted_outfiles, @splitted_stderrfiles );

  TEMP_FILE:
    foreach my $temp_file ( $splitted_files, $vcf_header_path, $vcf_body_path ) {
        gnu_rm(
            {
                filehandle  => $filehandle,
                infile_path => $temp_file,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                base_command         => $profile_base_command,
                dependency_method    => q{island_to_samples},
                case_id              => $case_id,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                job_id_chain         => q{MAIN},
                recipe_file_path     => $recipe_file_path,
                sample_ids_ref       => \@{ $active_parameter_href->{sample_ids} },
                submission_profile   => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _check_outfile_and_move {

## Function: Check if the vcf tag has already been created before moving the file in place
## Returns:
## Arguments: $filehandle        => Filehandle
##          : $outfile_path      => Outfile
##          : $temp_outfile_path => Temporary outfile
##          : $vcf_id_tag        => Recipe name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $outfile_path;
    my $reference_dir;
    my $temp_outfile_path;
    my $vcf_id_tag;

    my $tmpl = {
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        reference_dir => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir,
            strict_type => 1,
        },
        temp_outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$temp_outfile_path,
            strict_type => 1,
        },
        vcf_id_tag => {
            defined     => 1,
            required    => 1,
            store       => \$vcf_id_tag,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv gnu_rm };

    ## Open bash if statment
    print {$filehandle} $OPEN_BRACKET x 2 . $SPACE . $DOLLAR_SIGN . $OPEN_PARENTHESIS;

    ## Insert condition
    perl_nae_oneliners(
        {
            filehandle         => $filehandle,
            oneliner_name      => q{get_vcf_header_id_line},
            oneliner_parameter => $vcf_id_tag,
            stdinfile_path     => $outfile_path,
        }
    );

    ## Close bash if statment
    print {$filehandle} $CLOSE_PARENTHESIS . $SPACE . $CLOSE_BRACKET x 2 . $SPACE;

    ## In case of OK
    print {$filehandle} $AMPERSAND x 2 . $SPACE;

    ## Remove temporary file
    gnu_rm(
        {
            filehandle  => $filehandle,
            infile_path => $temp_outfile_path . $ASTERISK,
        }
    );

    ## Otherwise
    print {$filehandle} $PIPE x 2 . $SPACE;

    ## Move temporary file in place
    gnu_mv(
        {
            filehandle   => $filehandle,
            infile_path  => $temp_outfile_path . $ASTERISK,
            outfile_path => $reference_dir,
        }
    );
    return;
}

sub _build_file_splitted_file_paths {

## Function: Build the splitted file paths. Based on the number of cores
## Returns: @splitted_files
## Arguments: $file_path_prefix => File path prefix
##          : $file_suffix      => File suffix
##          : $max_file_number  => Max file number

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path_prefix;
    my $file_suffix;
    my $max_file_number;

    my $tmpl = {
        file_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$file_path_prefix,
            strict_type => 1,
        },
        file_suffix => {
            default     => $arg_href->{file_suffix} ||= $EMPTY_STR,
            store       => \$file_suffix,
            strict_type => 1,
        },
        max_file_number => {
            defined     => 1,
            required    => 1,
            store       => \$max_file_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @splitted_files = map { $file_path_prefix . $DOT . $_ . $file_suffix }
      ( q{00} .. "$max_file_number" );

    return @splitted_files;
}

1;
