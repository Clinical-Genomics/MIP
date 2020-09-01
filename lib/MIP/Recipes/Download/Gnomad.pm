package MIP::Recipes::Download::Gnomad;

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
use MIP::Constants qw{ $DASH $DOT $NEWLINE $PIPE $SINGLE_QUOTE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_gnomad };

}

sub download_gnomad {

## Function : Download gnomad
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
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_dead_end };
    use MIP::Program::Bcftools qw{ bcftools_annotate bcftools_index };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Script::Setup_script qw{ setup_script };

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

## Map of key names to keep from reference vcf
    my %info_key = (
        q{r2.0.1}    => [ qw{ INFO/AF INFO/AF_POPMAX }, ],
        q{r2.1.1}    => [ qw{ INFO/AF INFO/AF_popmax }, ],
        q{r2.1.1_sv} => [ qw{ INFO/AC INFO/AF INFO/POPMAX_AF }, ],
    );

    my $reformated_outfile = join $UNDERSCORE,
      (
        $genome_version, $recipe_name, q{reformated},
        $DASH . $reference_version . q{-.vcf.gz}
      );
    my $reformated_outfile_path = catfile( $reference_dir, $reformated_outfile );

    ## Only include sites for which at least one of the info keys are above zero
    my $include = join $SPACE . $PIPE x 2 . $SPACE,
      map { $_ . q{>0} } @{ $info_key{$reference_version} };
    bcftools_annotate(
        {
            filehandle     => $filehandle,
            include        => $SINGLE_QUOTE . $include . $SINGLE_QUOTE,
            infile_path    => catfile( $reference_dir, $reference_href->{outfile} ),
            outfile_path   => $reformated_outfile_path,
            output_type    => q{z},
            remove_ids_ref => [ map { q{^} . $_ } @{ $info_key{$reference_version} } ],
        }
    );
    say {$filehandle} $NEWLINE;

    bcftools_index(
        {
            filehandle  => $filehandle,
            infile_path => $reformated_outfile_path,
            output_type => q{tbi},
        }
    );
    say {$filehandle} $NEWLINE;

    ## Create AF file for bcftools roh
    _build_af_file(
        {
            filehandle        => $filehandle,
            file_name         => $reference_href->{outfile},
            infile_path       => $reformated_outfile_path,
            reference_dir     => $reference_dir,
            reference_version => $reference_version,
        }
    );

    ## Close filehandle
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

sub _build_af_file {

## Function : Build allele frequency file for bcftools roh
## Returns  :
## Arguments: $file_name         => Name of downloaded file
##          : $filehandle        => Filehandle
##          : $infile_path       => Path to reformatted file
##          : $reference_dir     => Reference directory
##          : $reference_version => Reference version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;
    my $filehandle;
    my $infile_path;
    my $reference_dir;
    my $reference_version;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        reference_dir => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir,
            strict_type => 1,
        },
        reference_version => {
            allow       => [qw{ r2.0.1 r2.1.1 r2.1.1_sv }],
            required    => 1,
            store       => \$reference_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parse::File qw{ parse_file_suffix };
    use MIP::Program::Bcftools qw{ bcftools_query };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };

    ## Don't build file for SV:s
    return if ( $reference_version eq q{r2.1.1_sv} );

    my $outfile_no_suffix = parse_file_suffix(
        {
            file_name   => $file_name,
            file_suffix => q{.vcf},
        }
    );
    my $allele_frq_file_path = catfile( $reference_dir, $outfile_no_suffix . q{.tab.gz} );

    bcftools_query(
        {
            filehandle       => $filehandle,
            format           => q{'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n'},
            infile_paths_ref => [$infile_path],
        }
    );
    print {$filehandle} $SPACE . $PIPE . $SPACE;

    htslib_bgzip(
        {
            filehandle      => $filehandle,
            stdoutfile_path => $allele_frq_file_path,
            write_to_stdout => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    htslib_tabix(
        {
            begin       => 2,
            end         => 2,
            filehandle  => $filehandle,
            infile_path => $allele_frq_file_path,
            sequence    => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

1;
