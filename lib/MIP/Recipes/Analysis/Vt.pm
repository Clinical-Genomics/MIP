package MIP::Recipes::Analysis::Vt;

use 5.026;
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
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $ASTERISK $DASH $DOT $EMPTY_STR $LOG_NAME $NEWLINE $PIPE $SEMICOLON $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vt analysis_vt_panel };

}

sub analysis_vt {

## Function : Split multi allelic records into single records and normalize
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_path               => File path
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_path;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Vt_core qw{ analysis_vt_core_rio};
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
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my %infile_path        = %{ $io{in}{file_path_href} };

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
    my $core_number = $recipe_resource{core_number};

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
                temp_directory   => $temp_directory,
            }
        )
    );

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
            core_number                     => $core_number,
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );
    my $stderr_path = $recipe_info_path . $DOT . q{stderr.txt};

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) = splitpath($stderr_path);

    ### SHELL:

    say {$filehandle}
q{## vt - Decompose (split multi allelic records into single records) and/or normalize variants};

    my $xargs_file_path_prefix;

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            filehandle         => $filehandle,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Split vcf into contigs
  CONTIG:
    while ( my ( $contig_index, $contig ) = each @contigs_size_ordered ) {

        ## vt - Split multi allelic records into single records and normalize
        analysis_vt_core_rio(
            {
                active_parameter_href  => $active_parameter_href,
                cmd_break              => $SEMICOLON,
                contig                 => $contig,
                decompose              => $active_parameter_href->{vt_decompose},
                filehandle             => $xargsfilehandle,
                gnu_sed                => 1,
                infile_path            => $infile_path{$contig},
                instream               => 0,
                normalize              => $active_parameter_href->{vt_normalize},
                outfile_path           => $outfile_path{$contig},
                uniq                   => $active_parameter_href->{vt_uniq},
                xargs_file_path_prefix => $xargs_file_path_prefix,
            }
        );

        if (   $contig_index == 0
            && $recipe_mode == 1 )
        {

            ## Split to enable submission to sample_info QC later
            my ( $volume_xargs, $directory_xargs, $stderr_file_xargs ) =
              splitpath($xargs_file_path_prefix);

            ## Collect QC metadata info for later use
            my $qc_vt_outfile_path = catfile( $directory,
                $stderr_file_xargs . $DOT . $contig . $DOT . q{stderr.txt} );
            set_recipe_outfile_in_sample_info(
                {
                    path             => $qc_vt_outfile_path,
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        ## VEP does not annotate '*' since the alt allele does not exist, this is captured in the upstream indel and SNV record associated with '*'
        ## Remove decomposed '*' entries
        if ( $active_parameter_href->{vt_missing_alt_allele} ) {

            # Update file tag
            my $alt_file_tag = $UNDERSCORE . q{nostar};

            my $removed_outfile_path =
              $outfile_path_prefix . $alt_file_tag . $DOT . $contig . $outfile_suffix;
            _remove_decomposed_asterisk_entries_xargs(
                {
                    contig                 => $contig,
                    infile_path            => $outfile_path{$contig},
                    outfile_path           => $removed_outfile_path,
                    xargsfilehandle        => $xargsfilehandle,
                    xargs_file_path_prefix => $xargs_file_path_prefix,
                }
            );

            gnu_mv(
                {
                    filehandle   => $xargsfilehandle,
                    infile_path  => $removed_outfile_path,
                    outfile_path => $outfile_path{$contig},
                }
            );
            say {$xargsfilehandle} $NEWLINE;

        }
    }

    say {$filehandle} q{wait}, $NEWLINE;
    close $filehandle or $log->logcroak(q{Could not close filehandle});
    close $xargsfilehandle
      or $log->logcroak(q{Could not close xargsfilehandle});

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                base_command            => $profile_base_command,
                case_id                 => $case_id,
                dependency_method       => q{sample_to_case},
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_chain            => $job_id_chain,
                job_id_href             => $job_id_href,
                job_reservation_name    => $active_parameter_href->{job_reservation_name},
                log                     => $log,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub analysis_vt_panel {

## Function : Split multi allelic records into single records and normalize
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_path               => File path
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_path;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Vt_core qw{ analysis_vt_core_rio};
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
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

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
                temp_directory         => $temp_directory,
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
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    say {$filehandle}
q{## vt - Decompose (split multi allelic records into single records) and/or normalize variants};

    ## vt - Split multi allelic records into single records and normalize
    analysis_vt_core_rio(
        {
            active_parameter_href => $active_parameter_href,
            cmd_break             => $SEMICOLON,
            decompose             => $active_parameter_href->{vt_decompose},
            filehandle            => $filehandle,
            gnu_sed               => 1,
            infile_path           => $infile_path,
            instream              => 0,
            normalize             => $active_parameter_href->{vt_normalize},
            outfile_path          => $outfile_path,
            uniq                  => $active_parameter_href->{vt_uniq},
        }
    );
    say {$filehandle} $NEWLINE;

    ## VEP does not annotate '*' since the alt allele does not exist, this is captured in the upstream indel and SNV record associated with '*'
    ## Remove decomposed '*' entries
    if ( $active_parameter_href->{vt_missing_alt_allele} ) {

        my $removed_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{nostar} . $outfile_suffix;
        _remove_decomposed_asterisk_entries(
            {
                filehandle   => $filehandle,
                infile_path  => $outfile_path,
                outfile_path => $removed_outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;

        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $removed_outfile_path,
                outfile_path => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    close $filehandle;

    if ( $recipe_mode == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                path             => $recipe_info_path . $DOT . q{stderr.txt},
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command            => $profile_base_command,
                case_id                 => $case_id,
                dependency_method       => q{sample_to_case},
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_chain            => $job_id_chain,
                job_id_href             => $job_id_href,
                job_reservation_name    => $active_parameter_href->{job_reservation_name},
                log                     => $log,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _remove_decomposed_asterisk_entries_xargs {

## Function : Remove decomposed '*' entries
## Returns  :
## Arguments: $contig                 => Contig
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $xargsfilehandle        => XARGS file handle
##          : $xargs_file_path_prefix => Xargs file path prefix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contig;
    my $infile_path;
    my $outfile_path;
    my $xargsfilehandle;
    my $xargs_file_path_prefix;

    my $tmpl = {
        contig => {
            defined  => 1,
            required => 1,
            store    => \$contig,
        },
        infile_path => {
            defined  => 1,
            required => 1,
            store    => \$infile_path,
        },
        outfile_path => {
            defined  => 1,
            required => 1,
            store    => \$outfile_path,
        },
        xargsfilehandle => {
            defined  => 1,
            required => 1,
            store    => \$xargsfilehandle,
        },
        xargs_file_path_prefix => {
            defined  => 1,
            required => 1,
            store    => \$xargs_file_path_prefix,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    my $remove_star_regexp = q?perl -nae \'?;

    ## Print if line does not contain asterisk
    $remove_star_regexp .= q?unless\($F\[4\] eq \"\*\") \{print $_\}\' ?;

    ## Print regexp
    print {$xargsfilehandle} $remove_star_regexp;

    ## Print infile
    print {$xargsfilehandle} $infile_path . $SPACE;

    ## Print outfile
    print {$xargsfilehandle} q{>} . $SPACE . $outfile_path . $SPACE;

    ## Print stderr file
    print {$xargsfilehandle} q{2>>}
      . $SPACE
      . $xargs_file_path_prefix
      . $DOT
      . $contig
      . $DOT
      . q{stderr.txt}
      . $SPACE;

    # Redirect xargs output to recipe specific stderr file
    print {$xargsfilehandle} $SEMICOLON . $SPACE;

    return;
}

sub _remove_decomposed_asterisk_entries {

## Function : Remove decomposed '*' entries
## Returns  :
## Arguments: $filehandle   => Filehandle
##          : $infile_path  => Infile path
##          : $outfile_path => Outfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;

    my $tmpl = {
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        infile_path => {
            defined  => 1,
            required => 1,
            store    => \$infile_path,
        },
        outfile_path => {
            defined  => 1,
            required => 1,
            store    => \$outfile_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    my $remove_star_regexp = q?perl -nae '?;

    ## Print if line does not contain asterisk
    $remove_star_regexp .= q?unless($F[4] eq q{*}) {print $_}' ?;

    ## Print regexp
    print {$filehandle} $remove_star_regexp;

    ## Print infile
    print {$filehandle} $infile_path . $SPACE;

    ## Print outfile
    print {$filehandle} q{>} . $SPACE . $outfile_path . $SPACE;

    return;
}

1;
