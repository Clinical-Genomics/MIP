package MIP::Recipes::Download::Cadd_whole_genome_snvs;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile catdir};
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
use MIP::Constants qw{ $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_cadd_whole_genome_snvs };

}

sub download_cadd_whole_genome_snvs {

## Function : Download cadd score for whole genome snvs
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
    use MIP::Gnu::Coreutils qw{ gnu_mkdir gnu_echo };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_dead_end };

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
    my $FILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href      => $active_parameter_href,
            core_number                => $recipe_resource{core_number},
            directory_id               => q{mip_download},
            FILEHANDLE                 => $FILEHANDLE,
            job_id_href                => $job_id_href,
            log                        => $log,
            memory_allocation          => $recipe_resource{memory},
            outdata_dir                => $reference_dir,
            outscript_dir              => $reference_dir,
            process_time               => $recipe_resource{time},
            recipe_data_directory_path => $active_parameter_href->{reference_dir},
            recipe_directory           => $recipe_name . $UNDERSCORE . $reference_version,
            recipe_name                => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### SHELL:

    say {$FILEHANDLE} q{## } . $recipe_name;

    ## Construct outdir path
    my $outdir_path =
      catdir( $reference_dir, qw{ CADD-scripts data prescored GRCh38_v1.5 no_anno } );
    if ( $genome_version eq q{grch37} and $reference_version eq q{v1.4} ) {

        $outdir_path =
          catdir( $reference_dir, qw{ CADD-scripts data prescored GRCh37_v1.4 no_anno } );
    }

    if ( not -d $outdir_path ) {
        gnu_mkdir(
            {
                FILEHANDLE       => $FILEHANDLE,
                indirectory_path => $outdir_path,
                parents          => 1,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    get_reference(
        {
            FILEHANDLE     => $FILEHANDLE,
            recipe_name    => $recipe_name,
            reference_dir  => $outdir_path,
            reference_href => $reference_href,
            quiet          => $quiet,
            verbose        => $verbose,
        }
    );

    ## Create some mock files for the later download processes to recognize
    my @file_paths = (
        catfile( $reference_dir, $reference_href->{outfile} ),
        catfile( $reference_dir, $reference_href->{outfile_check} ),
        catfile( $reference_dir, $reference_href->{outfile_index} ),
        catfile( $reference_dir, $reference_href->{outfile_index_check} ),
    );
    my $echo_message = q{File downloaded to} . $SPACE . $outdir_path;

  FILE_PATH:
    foreach my $file_path (@file_paths) {
        gnu_echo(
            {
                outfile_path => $file_path,
                FILEHANDLE   => $FILEHANDLE,
                strings_ref  => [$echo_message],
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }

    ## Close FILEHANDLES
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

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

1;
