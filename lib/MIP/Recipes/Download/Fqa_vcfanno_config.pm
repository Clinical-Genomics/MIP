package MIP::Recipes::Download::Fqa_vcfanno_config;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile };
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
  qw{ $DOT $DOUBLE_QUOTE $EMPTY_STR $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_fqa_vcfanno_config };

}

sub download_fqa_vcfanno_config {

## Function : Download fqa_vcfanno_config
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
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv};
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_dead_end };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
            active_parameter_href      => $active_parameter_href,
            core_number                => $recipe_resource{core_number},
            directory_id               => q{mip_download},
            filehandle                 => $filehandle,
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

    my $reference_dir_path = catfile( $reference_dir, $EMPTY_STR );
    my $outfile_path       = catfile( $reference_dir, $reference_href->{outfile} );
    my $outfile_path_tmp   = $outfile_path . $DOT . q{tmp};

    my $set_reference_path_regexp = _set_reference_dir_path(
        {
            infile_path        => $outfile_path,
            outfile_path       => $outfile_path_tmp,
            reference_dir_path => $reference_dir_path,
        }
    );

    say {$filehandle} $set_reference_path_regexp, $NEWLINE;

    gnu_mv(
        {
            filehandle   => $filehandle,
            infile_path  => $outfile_path_tmp,
            outfile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

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

sub _set_reference_dir_path {

## Function : Prepend the reference directory path to the reference file name in config
## Returns  : $set_reference_path_regexp
## Arguments: $infile_path        => Infile path
##          : $outfile_path       => Outfile path
##          : $reference_dir_path => Reference dir path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $reference_dir_path;

    my $tmpl = {
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
        reference_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    my $set_reference_path_regexp = q?perl -nae '?;

    ## Find "file" line  and catch file name between "
    $set_reference_path_regexp .= q? my ($filename) = $_=~/file=["](\S+)["]/sxm;?;

    ## If file name prepend reference dir to file name
    $set_reference_path_regexp .=
        q? if($filename) {say STDOUT q{?
      . q{file=}
      . $DOUBLE_QUOTE
      . $reference_dir_path
      . q?} . $filename . q{?
      . $DOUBLE_QUOTE . q?}}?;

    ## Else print other lines as is
    $set_reference_path_regexp .= q? else {print STDOUT $_} ' ?;

    ## Infile
    $set_reference_path_regexp .= $infile_path;

    ## Outfile
    $set_reference_path_regexp .= q? > ? . $outfile_path;

    return $set_reference_path_regexp;
}

1;
