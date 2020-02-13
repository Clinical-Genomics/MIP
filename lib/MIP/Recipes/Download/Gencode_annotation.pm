package MIP::Recipes::Download::Gencode_annotation;

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
use MIP::Constants qw{ $BACKWARD_SLASH $DASH $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_gencode_annotation };

}

sub download_gencode_annotation {

## Function : Download gencode_annotation
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
    use MIP::Program::Gtf2bed qw{ gtf2bed };
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

    ### POST PROCESSING
    my $outfile_name = join $UNDERSCORE,
      ( $genome_version, $recipe_name, q{-} . $reference_version . q{-.gtf} );
    my $outfile_path = catfile( $reference_dir, $outfile_name );

    if ( $genome_version eq q{grch37} ) {
        ## Build reformated outfile
        my $reformated_outfile = join $UNDERSCORE,
          (
            $genome_version, $recipe_name, q{reformated},
            q{-} . $reference_version . q{-.gtf}
          );
        my $reformated_outfile_path = catfile( $reference_dir, $reformated_outfile );

        _remove_chr_prefix_rename_chrm(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path,
            }
        );
        say {$filehandle} $PIPE . $SPACE . $BACKWARD_SLASH;

        _remove_chr_prefix(
            {
                filehandle   => $filehandle,
                outfile_path => $reformated_outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;

        $outfile_path = $reformated_outfile_path;
    }

    ## Reformat gtf to bed
    my $bed_outfile_path = $outfile_path =~ s/gtf$/bed/rxms;
    gtf2bed(
        {
            filehandle      => $filehandle,
            infile_path     => $outfile_path,
            stdoutfile_path => $bed_outfile_path,
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

sub _remove_chr_prefix_rename_chrm {

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

    # Execute perl
    print {$filehandle} q{perl -nae ' };

    # Print header
    print {$filehandle} q?if($_=~/^#/) { print $_;} ?;

    ## Remove prefix and rename mitochondria prefix
    print {$filehandle} q?else { $_ =~ s/^(chrM)/MT/g; print $_;}' ?;

    # Infile
    print {$filehandle} $infile_path . $SPACE;

    return;
}

sub _remove_chr_prefix {

## Function : Remove chr prefix in file
## Returns  :
## Arguments: $filehandle   => Filehandle to write to
##          : $outfile_path => Outfile path

    my ($arg_href) = @_;

## Flatten argument(s)
    my $filehandle;
    my $outfile_path;

    my $tmpl = {
        filehandle   => { defined => 1, required => 1, store => \$filehandle, },
        outfile_path => {
            default     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Remove chr prefix
    # Execute perl
    print {$filehandle} q{perl -nae ' };

    # Print header
    print {$filehandle} q?if($_=~/^#/) { print $_; } ?;

    # Remove prefix
    print {$filehandle} q?else {$_ =~ s/^chr(.+)/$1/g; print $_; }' ?;

    # Infile
    print {$filehandle} $DASH . $SPACE;

    # Outfile
    print {$filehandle} q{> } . $outfile_path;

    return;
}

1;
