package MIP::Recipes::Analysis::Vt_core;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
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
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.10;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vt_core analysis_vt_core_rio };
}

sub analysis_vt_core {

## Function : Split multi allelic records into single records and normalize
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $bcftools_output_type    => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $bgzip                   => Compress output from vt using bgzip
##          : $cmd_break               => Command line separator ['"\n\n"'|";"]
##          : $contig                  => The contig to extract {OPTIONAL, REF}
##          : $core_number             => The number of cores to allocate
##          : $decompose               => Vt program decompose for splitting multiallelic variants
##          : $case_id                 => The case ID
##          : $filehandle              => Filehandle to write to
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $infile_path             => Infile path
##          : $instream                => Data to vt is supplied as a unix pipe
##          : $gnu_sed                 => Sed program for changing vcf #FORMAT field in variant vcfs
##          : $human_genome_reference  => Human genome reference
##          : $normalize               => Vt program normalize for normalizing to reference used in analysis
##          : $outfile_path            => Outfile path
##          : $parameter_href          => Hash with paremters from yaml file {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_directory        => Program directory to write to in sbatch script
##          : $recipe_name             => Program name
##          : $tabix                   => Index compressed output using tabix
##          : $uniq                    => Vt program uniq for removing variant duplication that appear later in file
##          : $xargs_file_path_prefix  => The xargs sbatch script file name {OPTIONAL}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $contig;
    my $filehandle;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $infile_path;
    my $parameter_href;
    my $xargs_file_path_prefix;

    ## Default(s)
    my $bcftools_output_type;
    my $bgzip;
    my $cmd_break;
    my $core_number;
    my $decompose;
    my $case_id;
    my $gnu_sed;
    my $human_genome_reference;
    my $instream;
    my $normalize;
    my $outfile_path;
    my $profile_base_command;
    my $recipe_directory;
    my $recipe_name;
    my $tabix;
    my $uniq;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        bcftools_output_type => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$bcftools_output_type,
            strict_type => 1,
        },
        bgzip => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$bgzip,
            strict_type => 1,
        },
        cmd_break => { default => $NEWLINE x 2, store => \$cmd_break, strict_type => 1, },
        core_number => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        contig    => { store => \$contig, strict_type => 1, },
        decompose => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$decompose,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        gnu_sed    => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$gnu_sed,
            strict_type => 1,
        },
        human_genome_reference => {
            default     => $arg_href->{active_parameter_href}{human_genome_reference},
            store       => \$human_genome_reference,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        instream => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$instream,
            strict_type => 1,
        },
        job_id_href => { default => {}, store => \$job_id_href, strict_type => 1, },
        normalize   => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$normalize,
            strict_type => 1,
        },
        outfile_path => {
            default     => $arg_href->{infile_path},
            store       => \$outfile_path,
            strict_type => 1,
        },
        parameter_href => { default => {}, strict_type => 1, store => \$parameter_href, },
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_directory =>
          { default => q{vt_core}, store => \$recipe_directory, strict_type => 1, },
        recipe_name =>
          { default => q{vt_core}, store => \$recipe_name, strict_type => 1, },
        tabix => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$tabix,
            strict_type => 1,
        },
        uniq => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$uniq,
            strict_type => 1,
        },
        xargs_file_path_prefix =>
          { store => \$xargs_file_path_prefix, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Gnu::Software::Gnu_sed qw{ gnu_sed };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view bcftools_index };
    use MIP::Program::Vt qw{ vt_decompose vt_normalize vt_uniq };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Set MIP recipe name
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

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ## Create anonymous filehandle
    $filehandle = IO::Handle->new();

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
            recipe_name                     => $recipe_name,
            recipe_directory                => $recipe_directory,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### Split multi allelic records into single records and normalize

    ## Get parameters
    my $stderrfile_path;

    ## Write stderr for xargs process
    if ( defined $xargs_file_path_prefix && defined $contig ) {

        ## Redirect xargs output to program specific stderr file
        $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
    }

    ## Initate processing
    if ( not $instream ) {

        bcftools_view(
            {
                filehandle  => $filehandle,
                infile_path => $infile_path,
                output_type => $bcftools_output_type,
            }
        );
    }
    ## Replace #FORMAT field prior to smart decomposition (variant vcfs)
    if ($gnu_sed) {

        ## Pipe
        print {$filehandle} $PIPE . $SPACE;

        gnu_sed(
            {
                filehandle => $filehandle,
                script     => q{'s/ID=AD,Number=./ID=AD,Number=R/'},
            }
        );
    }
    if ($decompose) {

        ## Pipe
        print {$filehandle} $PIPE . $SPACE;

        vt_decompose(
            {
                filehandle             => $filehandle,
                infile_path            => q{-},
                smart_decomposition    => 1,
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path,
            }
        );
    }
    if ($normalize) {

        # Pipe
        print {$filehandle} $PIPE . $SPACE;

        vt_normalize(
            {
                filehandle                     => $filehandle,
                infile_path                    => q{-},
                no_fail_inconsistent_reference => 1,
                referencefile_path             => $human_genome_reference,
                stderrfile_path                => $stderrfile_path,
                stderrfile_path_append         => $stderrfile_path,
            }
        );
    }
    if ($uniq) {

        ## Pipe
        print {$filehandle} $PIPE . $SPACE;

        vt_uniq(
            {
                filehandle             => $filehandle,
                infile_path            => q{-},
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path,
            }
        );
    }

    ## Tabix has been/will be used on file, compress again
    my $tabix_suffix    = $DOT . q{tbi};
    my $infile_path_tbi = $infile_path . $tabix_suffix;
    if ( -e $infile_path_tbi || $bgzip ) {

        ## Pipe
        print {$filehandle} $PIPE . $SPACE;

        htslib_bgzip(
            {
                filehandle      => $filehandle,
                write_to_stdout => 1,
            }
        );
    }

    ## Temporary outfile
    print {$filehandle} q{>}
      . $SPACE
      . $outfile_path
      . $UNDERSCORE
      . q{splitted}
      . $UNDERSCORE
      . $random_integer
      . $SPACE;
    print {$filehandle} $cmd_break;

    ## tabix index
    if ( -e $infile_path_tbi || $tabix ) {

        my $tabix_infile_path =
          $outfile_path . $UNDERSCORE . q{splitted} . $UNDERSCORE . $random_integer;
        htslib_tabix(
            {
                filehandle  => $filehandle,
                force       => 1,
                infile_path => $tabix_infile_path,
                preset      => q{vcf},
            }
        );
        print {$filehandle} $cmd_break;

        ## Move index in place
        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $tabix_infile_path . $tabix_suffix,
                outfile_path => $outfile_path . $tabix_suffix,
            }
        );
        print {$filehandle} $cmd_break;
    }

    ## Move processed reference to original place
    gnu_mv(
        {
            filehandle  => $filehandle,
            infile_path => $outfile_path
              . $UNDERSCORE
              . q{splitted}
              . $UNDERSCORE
              . $random_integer,
            outfile_path => $outfile_path,
        }
    );
    print {$filehandle} $cmd_break;

    close $filehandle;

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

sub analysis_vt_core_rio {

## Function : Split multi allelic records into single records and normalize
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $bcftools_output_type    => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $bgzip                   => Compress output from vt using bgzip
##          : $cmd_break               => Command line separator ['"\n\n"'|";"]
##          : $contig                  => The contig to extract {OPTIONAL, REF}
##          : $core_number             => The number of cores to allocate
##          : $decompose               => Vt program decompose for splitting multiallelic variants
##          : $case_id                 => The case ID
##          : $filehandle              => Filehandle to write to
##          : $gnu_sed                 => Sed program for changing vcf #FORMAT field in variant vcfs
##          : $human_genome_reference  => Human genome reference
##          : $infile_path             => Infile path
##          : $instream                => Data to vt is supplied as a unix pipe
##          : $normalize               => Vt program normalize for normalizing to reference used in analysis
##          : $outfile_path            => Outfile path
##          : $recipe_directory        => Program directory to write to in sbatch script
##          : $recipe_name             => Program name
##          : $tabix                   => Index compressed output using tabix
##          : $uniq                    => Vt program uniq for removing variant duplication that appear later in file
##          : $xargs_file_path_prefix  => The xargs sbatch script file name {OPTIONAL}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $contig;
    my $filehandle;
    my $infile_path;
    my $xargs_file_path_prefix;

    ## Default(s)
    my $bcftools_output_type;
    my $bgzip;
    my $cmd_break;
    my $core_number;
    my $decompose;
    my $case_id;
    my $gnu_sed;
    my $human_genome_reference;
    my $instream;
    my $normalize;
    my $outfile_path;
    my $recipe_directory;
    my $recipe_name;
    my $tabix;
    my $uniq;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        bcftools_output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$bcftools_output_type,
        },
        bgzip => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$bgzip
        },
        cmd_break => { default => $NEWLINE x 2, strict_type => 1, store => \$cmd_break },
        contig      => { strict_type => 1, store => \$contig },
        core_number => {
            default     => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$core_number
        },
        decompose => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$decompose
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            strict_type => 1,
            store       => \$case_id
        },
        filehandle => { required => 1, store => \$filehandle },
        gnu_sed    => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$gnu_sed
        },
        human_genome_reference => {
            default     => $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        instream => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$instream
        },
        normalize => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$normalize
        },
        ## Use same path as infile path unless parameter is supplied
        outfile_path => {
            default     => $arg_href->{infile_path},
            strict_type => 1,
            store       => \$outfile_path,
        },
        recipe_directory =>
          { default => q{vt}, strict_type => 1, store => \$recipe_directory },
        recipe_name => { default => q{vt}, strict_type => 1, store => \$recipe_name },
        xargs_file_path_prefix => { strict_type => 1, store => \$xargs_file_path_prefix },
        tabix                  => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$tabix
        },
        uniq => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$uniq
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Gnu::Software::Gnu_sed qw{ gnu_sed };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_add_to_samples };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view bcftools_index };
    use MIP::Program::Vt qw{ vt_decompose vt_normalize vt_uniq };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Set MIP recipe name
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    my $recipe_info_path;

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ### Split multi allelic records into single records and normalize

    ## Get parameters
    my $stderrfile_path;

    ## Write stderr for xargs process
    if ( defined $xargs_file_path_prefix && defined $contig ) {

        ## Redirect xargs output to program specific stderr file
        $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
    }

    ## Initate processing
    if ( not $instream ) {

        bcftools_view(
            {
                filehandle  => $filehandle,
                infile_path => $infile_path,
                output_type => $bcftools_output_type,
            }
        );
    }
    ## Replace #FORMAT field prior to smart decomposition (variant vcfs)
    if ($gnu_sed) {

        ## Pipe
        print {$filehandle} $PIPE . $SPACE;

        gnu_sed(
            {
                filehandle => $filehandle,
                script     => q{'s/ID=AD,Number=./ID=AD,Number=R/'},
            }
        );
    }
    if ($decompose) {

        ## Pipe
        print {$filehandle} $PIPE . $SPACE;

        vt_decompose(
            {
                filehandle             => $filehandle,
                infile_path            => q{-},
                smart_decomposition    => 1,
                stderrfile_path_append => $stderrfile_path,
            }
        );
    }
    if ($normalize) {

        # Pipe
        print {$filehandle} $PIPE . $SPACE;

        vt_normalize(
            {
                filehandle                     => $filehandle,
                infile_path                    => q{-},
                no_fail_inconsistent_reference => 1,
                referencefile_path             => $human_genome_reference,
                stderrfile_path_append         => $stderrfile_path,
            }
        );
    }
    if ( $active_parameter_href->{vt_uniq} > 0 && $uniq ) {

        ## Pipe
        print {$filehandle} $PIPE . $SPACE;

        vt_uniq(
            {
                filehandle             => $filehandle,
                infile_path            => q{-},
                stderrfile_path_append => $stderrfile_path,
            }
        );
    }

    ## tabix has been/will be used on file, compress again
    my $tabix_suffix    = $DOT . q{tbi};
    my $infile_path_tbi = $infile_path . $tabix_suffix;
    if ( -e $infile_path_tbi || $bgzip ) {

        ## Pipe
        print {$filehandle} $PIPE . $SPACE;

        htslib_bgzip(
            {
                filehandle      => $filehandle,
                write_to_stdout => 1,
            }
        );
    }

    ## Temporary outfile
    print {$filehandle} q{>}
      . $SPACE
      . $outfile_path
      . $UNDERSCORE
      . q{splitted}
      . $UNDERSCORE
      . $random_integer
      . $SPACE;
    print {$filehandle} $cmd_break;

    ## tabix index
    if ( -e $infile_path_tbi || $tabix ) {

        my $tabix_infile_path =
          $outfile_path . $UNDERSCORE . q{splitted} . $UNDERSCORE . $random_integer;
        htslib_tabix(
            {
                filehandle  => $filehandle,
                force       => 1,
                infile_path => $tabix_infile_path,
                preset      => q{vcf},
            }
        );
        print {$filehandle} $cmd_break;

        ## Move index in place
        gnu_mv(
            {
                filehandle   => $filehandle,
                infile_path  => $tabix_infile_path . $tabix_suffix,
                outfile_path => $outfile_path . $tabix_suffix,
            }
        );
        print {$filehandle} $cmd_break;
    }

    ## Move processed reference to original place
    gnu_mv(
        {
            filehandle  => $filehandle,
            infile_path => $outfile_path
              . $UNDERSCORE
              . q{splitted}
              . $UNDERSCORE
              . $random_integer,
            outfile_path => $outfile_path,
        }
    );
    print {$filehandle} $cmd_break;

    return 1;
}

1;
