package MIP::Recipes::Download::Clinvar;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ remove_tree };
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
use MIP::Constants qw{ $BACKWARD_SLASH $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_clinvar };

}

sub download_clinvar {

## Function : Download clinvar references
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
    use MIP::Gnu::Coreutils qw{ gnu_rm };
    use MIP::Program::Utility::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_annotate };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_dead_end };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Unpack parameters
    my $reference_dir = $active_parameter_href->{reference_dir};
    my @reference_genome_versions =
      @{ $active_parameter_href->{reference_genome_versions} };
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
            temp_directory             => $temp_directory,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### SHELL:

    say {$FILEHANDLE} q{## } . $recipe_name;

    get_reference(
        {
            FILEHANDLE     => $FILEHANDLE,
            recipe_name    => $recipe_name,
            reference_dir  => $reference_dir,
            reference_href => $reference_href,
            quiet          => $quiet,
            verbose        => $verbose,
        }
    );

    say {$FILEHANDLE} q{## Build clinvar variation ID header file};
    my $header_file_path =
      catfile( $reference_dir, $reference_version . $UNDERSCORE . q{clnvid_header.txt} );
    ## Build clinvar variation ID header file
    _build_clnvid_head_file(
        {
            FILEHANDLE       => $FILEHANDLE,
            header_file_path => $header_file_path,
        }
    );

    say {$FILEHANDLE} q{## Annotate vcf with new header};
    bcftools_annotate(
        {
            FILEHANDLE      => $FILEHANDLE,
            headerfile_path => $header_file_path,
            infile_path     => catfile( $reference_dir, $reference_href->{outfile} ),
            output_type     => q{v},
        }
    );
    say {$FILEHANDLE} $PIPE . $SPACE . $BACKWARD_SLASH;

    my $reformated_outfile = join $UNDERSCORE,
      (
        $genome_version, $recipe_name, q{reformated}, q{-} . $reference_version . q{-.vcf}
      );
    my $reformated_outfile_path = catfile( $reference_dir, $reformated_outfile );

    ## Add clinvar variation ID to vcf info file
    _add_clnvid_to_vcf_info(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $reformated_outfile_path,
        }
    );

    say {$FILEHANDLE} q{## Compress and index file};
    ## Compress file
    htslib_bgzip(
        {
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
            infile_path => $reformated_outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Index file using tabix
    htslib_tabix(
        {
            FILEHANDLE  => $FILEHANDLE,
            force       => 1,
            infile_path => $reformated_outfile_path . q{.gz},
            preset      => q{vcf},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Clean-up
    gnu_rm(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $header_file_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

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

sub _build_clnvid_head_file {

    ## Function : Build clinvar variation ID header file
    ## Returns  :
    ## Arguments: $FILEHANDLE       => Filehandle to write to
    ##          : $header_file_path => VCF header file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $header_file_path;

    my $tmpl = {
        FILEHANDLE       => { store => \$FILEHANDLE, },
        header_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$header_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    print {$FILEHANDLE} q{perl -e '};

    ## Print header line for Clinvar variation ID
    print {$FILEHANDLE}
q? print q{##INFO=<ID=CLNVID,Number=1,Type=Integer,Description="ClinVar Variation ID">} '?;

    ## Write to files
    say {$FILEHANDLE} q{ > } . $header_file_path . $NEWLINE;

    return;
}

sub _add_clnvid_to_vcf_info {

    ## Function : Add clinvar variation ID to vcf info file
    ## Returns  :
    ## Arguments: $FILEHANDLE   => Filehandle to write to
    ##          : $outfile_path => VCF header file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $outfile_path;

    my $tmpl = {
        FILEHANDLE   => { store => \$FILEHANDLE, },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    print {$FILEHANDLE} q{perl -nae ' };

    ## Skip header lines
    print {$FILEHANDLE} q?if($_=~/^#/) { print $_;} ?;

    ## Else add CLVID to INFO
    print {$FILEHANDLE}
      q?else { chomp; my $line = $_; say STDOUT $_ . q{;CLNVID=} . $F[2] } '?;

    ## Write to files
    say {$FILEHANDLE} q{ > } . $outfile_path . $NEWLINE;

    return;
}

1;
