package MIP::Recipes::Download::Dbnsfp;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $BACKWARD_SLASH $COMMA $DASH $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_dbnsfp };

}

## Constants
Readonly my $GRCH38_CHR_POS    => 1;
Readonly my $GRCH38_REGION_POS => 2;
Readonly my $GRCH37_CHR_POS    => 8;
Readonly my $GRCH37_REGION_POS => 9;

sub download_dbnsfp {

## Function : Download dbnsfp
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $genome_version        => Human genome version
##          : $job_id_href           => The job_id hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Recipe name
##          : $reference_href        => Reference hash {REF}
##          : $reference_version     => Reference version
##          : $quiet                 => Quiet (no output)
##          : $temp_directory        => Temporary directory for recipe
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $genome_version;
    my $job_id_href;
    my $recipe_name;
    my $reference_href;
    my $reference_version;

    ## Default(s)
    my $profile_base_command;
    my $quiet;
    my $temp_directory;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        genome_version => {
            store       => \$genome_version,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
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
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
        reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cat gnu_head gnu_sort };
    use MIP::Program::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::Language::Awk qw{ awk };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Processmanagement::Slurm_processes qw{ slurm_submit_job_no_dependency_dead_end };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Unpack parameters
    my $reference_dir = $active_parameter_href->{reference_dir};

    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set recipe mode
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    ## Filehandle(s)
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => q{mip_download},
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            outdata_dir                     => $reference_dir,
            outscript_dir                   => $reference_dir,
            process_time                    => $recipe_resource{time},
            recipe_data_directory_path      => $active_parameter_href->{reference_dir},
            recipe_directory                => $recipe_name . $UNDERSCORE . $reference_version,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    get_reference(
        {
            filehandle     => $filehandle,
            recipe_name    => $recipe_name,
            reference_dir  => $reference_dir,
            reference_href => $reference_href,
            quiet          => $quiet,
            verbose        => $verbose,
        }
    );

    say {$filehandle} q{## Build dbnsfp header file after unzip};

    ## Build dbnsfp chr file name after unzip
    my $dbnsfp_chr_file_name = q{dbNSFP} . $reference_version . $UNDERSCORE . q{variant.chr};
    my $dbnsfp_chr_file_path = catfile( $reference_dir, $dbnsfp_chr_file_name );
    my $reformated_outfile   = join $UNDERSCORE,
      ( $genome_version, $recipe_name, q{reformated}, $DASH . $reference_version . q{-.txt} );
    my $reformated_outfile_path = catfile( $reference_dir, $reformated_outfile );

    if ( $genome_version eq q{grch37} ) {

        ## Switch columns to use with genome build prior to version 38
        _reformat_for_grch37(
            {
                filehandle        => $filehandle,
                infile_path       => $dbnsfp_chr_file_path,
                outfile_path      => $reformated_outfile_path,
                reference_dir     => $reference_dir,
                reference_version => $reference_version,
                temp_directory    => $temp_directory,
            }
        );
    }
    else {

## Reformat by bgzip and tabix index
        _reformat_dbnsfp(
            {
                filehandle        => $filehandle,
                infile_path       => $dbnsfp_chr_file_path,
                outfile_path      => $reformated_outfile_path,
                reference_dir     => $reference_dir,
                reference_version => $reference_version,
                temp_directory    => $temp_directory,
            }
        );
    }

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        ## No upstream or downstream dependencies
        slurm_submit_job_no_dependency_dead_end(
            {
                base_command     => $profile_base_command,
                job_id_href      => $job_id_href,
                log              => $log,
                sbatch_file_name => $recipe_file_path,
            }
        );
    }
    return 1;
}

sub _reformat_for_grch37 {

## Function : Switch columns to use with genome build prior to version 38
## Returns  :
## Arguments: $filehandle        => Filehandle to write to
##          : $infile_path       => Infile path
##          : $outfile_path      => Outfile path
##          : $reference_dir     => Reference directory
##          : $reference_version => Reference version
##          : $temp_directory    => Temporary directory for recipe

    my ($arg_href) = @_;

## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $reference_dir;
    my $reference_version;
    my $temp_directory;

    my $tmpl = {
        filehandle  => { defined => 1, required => 1, store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
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
        reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $dbnsfp_chr1_file_path = $infile_path . q{1};

    my $header_file_path = catfile( $reference_dir, q{dbnsfp_header.txt} );

    say {$filehandle} q{## Switch columns to use with genome build prior to version 38};
    gnu_head(
        {
            filehandle      => $filehandle,
            infile_path     => $dbnsfp_chr1_file_path,
            lines           => 1,
            stdoutfile_path => $header_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Read files
    gnu_cat(
        {
            filehandle       => $filehandle,
            infile_paths_ref => [ $infile_path . $ASTERISK ],
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Skip header line starting with "#chr"
    gnu_grep(
        {
            filehandle   => $filehandle,
            infile_path  => $DASH,
            invert_match => 1,
            pattern      => q{^#chr},
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## If not DOT in column 8
    awk(
        {
            filehandle => $filehandle,
            statement  => q{$8 != "."},
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Sort nummerical on column number for $GRCH37_CHR_POS and $GRCH37_REGION_POS
    my $sort_chr_pos    = join $COMMA, ($GRCH37_CHR_POS) x 2;
    my $sort_region_pos = join( $COMMA, ($GRCH37_REGION_POS) x 2 ) . q{n};
    gnu_sort(
        {
            filehandle          => $filehandle,
            infile_path         => $DASH,
            keys_ref            => [ $sort_chr_pos, $sort_region_pos ],
            temporary_directory => $temp_directory,
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Concatenate header and input stream (DASH)
    gnu_cat(
        {
            filehandle       => $filehandle,
            infile_paths_ref => [ $header_file_path, $DASH ],
            stdoutfile_path  => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Compress and index file};

    ## Compress file
    htslib_bgzip(
        {
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Index file using tabix
    htslib_tabix(
        {
            begin       => $GRCH37_REGION_POS,
            end         => $GRCH37_REGION_POS,
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $outfile_path . q{.gz},
            sequence    => $GRCH37_CHR_POS,
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

sub _reformat_dbnsfp {

## Function : Reformat by bgzip and tabix index
## Returns  :
## Arguments: $filehandle        => Filehandle to write to
##          : $infile_path       => Infile path
##          : $outfile_path      => Outfile path
##          : $reference_dir     => Reference directory
##          : $reference_version => Reference version
##          : $temp_directory    => Temporary directory for recipe

    my ($arg_href) = @_;

## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $reference_dir;
    my $reference_version;
    my $temp_directory;

    my $tmpl = {
        filehandle  => { defined => 1, required => 1, store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
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
        reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gzip qw{ gzip };

    ## Any chr will do to generate header info
    my $dbnsfp_chr1_file_path = $infile_path . q{1.gz};

    ## Header file
    my $header_file_path = catfile( $reference_dir, q{dbnsfp_header.txt} );

    ## Read chr file
    gzip(
        {
            decompress       => 1,
            filehandle       => $filehandle,
            infile_paths_ref => [$dbnsfp_chr1_file_path],
            stdout           => 1,
        }
    );
    print {$filehandle} $PIPE . $SPACE;

    gnu_grep(
        {
            filehandle      => $filehandle,
            pattern         => q{^#chr},
            infile_path     => $DASH,
            stdoutfile_path => $header_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Read files
    gzip(
        {
            decompress       => 1,
            filehandle       => $filehandle,
            force            => 1,
            infile_paths_ref => [ $infile_path . $ASTERISK ],
            stdout           => 1,
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Skip header line starting with "#chr"
    gnu_grep(
        {
            filehandle   => $filehandle,
            infile_path  => $DASH,
            invert_match => 1,
            pattern      => q{^#chr},
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    my $sort_chr_pos    = join $COMMA, ($GRCH38_CHR_POS) x 2;
    my $sort_region_pos = join( $COMMA, ($GRCH38_REGION_POS) x 2 ) . q{n};
    gnu_sort(
        {
            filehandle          => $filehandle,
            infile_path         => $DASH,
            keys_ref            => [ $sort_chr_pos, $sort_region_pos ],
            temporary_directory => $temp_directory,
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Concatenate header and input stream (DASH)
    gnu_cat(
        {
            filehandle       => $filehandle,
            infile_paths_ref => [ $header_file_path, $DASH ],
            stdoutfile_path  => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Compress and index file};

    ## Compress file
    htslib_bgzip(
        {
            filehandle  => $filehandle,
            infile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Index file using tabix
    htslib_tabix(
        {
            begin       => $GRCH38_REGION_POS,
            end         => $GRCH38_REGION_POS,
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $outfile_path . q{.gz},
            sequence    => $GRCH38_CHR_POS,
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}
1;
