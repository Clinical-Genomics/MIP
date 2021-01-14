package MIP::Recipes::Download::Genomic_superdups;

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
use MIP::Constants
  qw{ $BACKWARD_SLASH $DASH $DOT $NEWLINE $PIPE $SEMICOLON $SPACE $TAB $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_genomic_superdups };

}

## Constants
Readonly my $GENOME_BUILD_VERSION_37 => q{grch37};

sub download_genomic_superdups {

## Function : Download genomic_superdups
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

    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cut gnu_sort gnu_uniq };
    use MIP::Program::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::File::Path qw{ remove_file_path_suffix };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Processmanagement::Slurm_processes qw{ slurm_submit_job_no_dependency_dead_end };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Unpack parameters
    my $reference_dir = $active_parameter_href->{reference_dir};

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

## Filehandle(s)
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe{core_number},
            directory_id                    => q{mip_download},
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe{memory},
            outdata_dir                     => $reference_dir,
            outscript_dir                   => $reference_dir,
            process_time                    => $recipe{time},
            recipe_data_directory_path      => $active_parameter_href->{reference_dir},
            recipe_directory                => $recipe_name . $UNDERSCORE . $reference_version,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe{load_env_ref},
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

    ## Parse file suffix in filename.suffix(.gz).
    ## Removes suffix if matching else return undef
    my $outfile_path_no_suffix = remove_file_path_suffix(
        {
            file_path         => catfile( $reference_dir, $reference_href->{outfile} ),
            file_suffixes_ref => [ $DOT . q{gz} ],
        }
    );

    ## Build reformated outfile
    my $reformated_outfile = join $UNDERSCORE,
      ( $genome_version, $recipe_name, q{reformated}, q{-} . $reference_version . q{-.bed} );
    my $reformated_outfile_path = catfile( $reference_dir, $reformated_outfile );

    if ( $genome_version eq $GENOME_BUILD_VERSION_37 ) {

        ## Repeat inline replace due to three chr entries  in file and use if "$1"
        say {$filehandle}
          q{## Repeat inline replace due to three chr entries  in file and use if "$1"};
        say   {$filehandle} q?for nr_chr_in_line in {1..3}?;
        say   {$filehandle} q{do};
        print {$filehandle} $TAB;
## Remove chr prefix in file
        _remove_chr_prefix(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path_no_suffix,
            }
        );
        say {$filehandle} $SEMICOLON;
        say {$filehandle} q{done} . $NEWLINE;

    }

    ## Reformat to bed
    say {$filehandle} q{## Reformat to bed};
    gnu_grep(
        {
            filehandle   => $filehandle,
            infile_path  => $outfile_path_no_suffix,
            invert_match => 1,
            pattern      => q{^#chr},
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Sort nummerical on column 2 and 3
    gnu_sort(
        {
            filehandle  => $filehandle,
            infile_path => $DASH,
            keys_ref    => [ q{2,2}, q{3,3n} ],
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Keep chr, start, stop and fraction match
    gnu_cut(
        {
            filehandle  => $filehandle,
            infile_path => $DASH,
            list        => q{2-4,27},
        }
    );
    say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

    ## Keep uniq entries only
    gnu_uniq(
        {
            filehandle      => $filehandle,
            infile_path     => $DASH,
            stdoutfile_path => $reformated_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Compress and index file};
    ## Compress file
    htslib_bgzip(
        {
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $reformated_outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Index file using tabix
    htslib_tabix(
        {
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $reformated_outfile_path . q{.gz},
            preset      => q{bed},
            zero_based  => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

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

sub _remove_chr_prefix {

## Function : Remove chr prefix in file
## Returns  :
## Arguments: $filehandle  => Filehandle to write to
##          : $infile_path => Infile path

    my ($arg_href) = @_;

## Flatten argument(s)
    my $filehandle;
    my $infile_path;

    my $tmpl = {
        filehandle  => { defined => 1, required => 1, store => \$filehandle, },
        infile_path => {
            default     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

## Remove chr prefix
    # Execute perl and modify file inline
    print {$filehandle} q{perl -i -p -e ' };

    # Unless header - remove prefix
    print {$filehandle} q?if($_!~/^#/) {s/chr(.+)/$1/g}' ?;

    # Infile to modify
    print {$filehandle} $infile_path;

    return;
}

1;
