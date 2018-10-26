package MIP::Recipes::Analysis::Vt_core;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vt_core analysis_vt_core_rio };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

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
##          : $case_id               => The case ID
##          : $FILEHANDLE              => Filehandle to write to
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $infile_path             => Infile path
##          : $instream                => Data to vt is supplied as a unix pipe
##          : $gnu_sed                 => Sed program for changing vcf #FORMAT field in variant vcfs
##          : $human_genome_reference  => Human genome reference
##          : $normalize               => Vt program normalize for normalizing to reference used in analysis
##          : $outfile_path            => Outfile path
##          : $parameter_href          => Hash with paremters from yaml file {REF}
##          : $recipe_directory       => Program directory to write to in sbatch script
##          : $recipe_name            => Program name
##          : $tabix                   => Index compressed output using tabix
##          : $uniq                    => Vt program uniq for removing variant duplication that appear later in file
##          : $xargs_file_path_prefix  => The xargs sbatch script file name {OPTIONAL}

    my ($arg_href) = @_;

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

    ## Flatten argument(s)
    my $active_parameter_href;
    my $contig;
    my $FILEHANDLE;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $infile_path;
    my $parameter_href;
    my $xargs_file_path_prefix;

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
        core_number => {
            default     => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$core_number
        },
        contig    => { strict_type => 1, store => \$contig },
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
        FILEHANDLE => { store => \$FILEHANDLE },
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
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
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
        job_id_href => { default => {}, strict_type => 1, store => \$job_id_href },
        ## Use same path as infile path unless parameter is supplied
        normalize => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$normalize
        },
        outfile_path => {
            default     => $arg_href->{infile_path},
            strict_type => 1,
            store       => \$outfile_path,
        },
        parameter_href => { default => {}, strict_type => 1, store => \$parameter_href },
        recipe_directory =>
          { default => q{vt}, strict_type => 1, store => \$recipe_directory },
        recipe_name => { default => q{vt}, strict_type => 1, store => \$recipe_name },
        tabix       => {
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
        xargs_file_path_prefix => { strict_type => 1, store => \$xargs_file_path_prefix },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Gnu::Software::Gnu_sed qw{ gnu_sed };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_add_to_samples };
    use MIP::Program::Utility::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Variantcalling::Allele_frequency qw{ calculate_af max_af };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view bcftools_index };
    use MIP::Program::Variantcalling::Vt qw{ vt_decompose vt_normalize vt_uniq };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;
    Readonly my $PROCESS_TIME      => 20;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP recipe name
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$recipe_name}{chain};

    my $file_path;
    my $recipe_info_path;

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ## Create anonymous filehandle
    $FILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    ( $file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $case_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            log                   => $log,
            process_time          => $PROCESS_TIME,
            recipe_name           => $recipe_name,
            recipe_directory      => $recipe_directory,
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
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $infile_path,
                output_type => $bcftools_output_type,
            }
        );
    }
    ## Replace #FORMAT field prior to smart decomposition (variant vcfs)
    if ($gnu_sed) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        gnu_sed(
            {
                FILEHANDLE => $FILEHANDLE,
                script     => q{'s/ID=AD,Number=./ID=AD,Number=R/'},
            }
        );
    }
    if ($decompose) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        vt_decompose(
            {
                FILEHANDLE             => $FILEHANDLE,
                infile_path            => q{-},
                smart_decomposition    => 1,
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path,
            }
        );
    }
    if ($normalize) {

        # Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        vt_normalize(
            {
                FILEHANDLE                     => $FILEHANDLE,
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
        print {$FILEHANDLE} $PIPE . $SPACE;

        vt_uniq(
            {
                FILEHANDLE             => $FILEHANDLE,
                infile_path            => q{-},
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path,
            }
        );
    }

    ## tabix has been/will be used on file, compress again
    my $tabix_suffix    = $DOT . q{tbi};
    my $infile_path_tbi = $infile_path . $tabix_suffix;
    if ( -e $infile_path_tbi || $bgzip ) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        htslib_bgzip(
            {
                FILEHANDLE      => $FILEHANDLE,
                write_to_stdout => 1,
            }
        );
    }

    ## Temporary outfile
    print {$FILEHANDLE} q{>}
      . $SPACE
      . $outfile_path
      . $UNDERSCORE
      . q{splitted}
      . $UNDERSCORE
      . $random_integer
      . $SPACE;
    print {$FILEHANDLE} $cmd_break;

    ## tabix index
    if ( -e $infile_path_tbi || $tabix ) {

        my $tabix_infile_path =
          $outfile_path . $UNDERSCORE . q{splitted} . $UNDERSCORE . $random_integer;
        htslib_tabix(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => $tabix_infile_path,
                preset      => q{vcf},
            }
        );
        print {$FILEHANDLE} $cmd_break;

        ## Move index in place
        gnu_mv(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $tabix_infile_path . $tabix_suffix,
                outfile_path => $outfile_path . $tabix_suffix,
            }
        );
        print {$FILEHANDLE} $cmd_break;
    }

    ## Move processed reference to original place
    gnu_mv(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path
              . $UNDERSCORE
              . q{splitted}
              . $UNDERSCORE
              . $random_integer,
            outfile_path => $outfile_path,
        }
    );
    print {$FILEHANDLE} $cmd_break;

    close $FILEHANDLE;

    if ( $recipe_mode == 1 ) {

        slurm_submit_job_no_dependency_add_to_samples(
            {
                case_id          => $case_id,
                job_id_href      => $job_id_href,
                log              => $log,
                path             => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
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
##          : $case_id               => The case ID
##          : $FILEHANDLE              => Filehandle to write to
##          : $gnu_sed                 => Sed program for changing vcf #FORMAT field in variant vcfs
##          : $human_genome_reference  => Human genome reference
##          : $infile_path             => Infile path
##          : $instream                => Data to vt is supplied as a unix pipe
##          : $normalize               => Vt program normalize for normalizing to reference used in analysis
##          : $outfile_path            => Outfile path
##          : $recipe_directory       => Program directory to write to in sbatch script
##          : $recipe_name            => Program name
##          : $tabix                   => Index compressed output using tabix
##          : $uniq                    => Vt program uniq for removing variant duplication that appear later in file
##          : $xargs_file_path_prefix  => The xargs sbatch script file name {OPTIONAL}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $contig;
    my $FILEHANDLE;
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
        FILEHANDLE => { store => \$FILEHANDLE },
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
    use MIP::Program::Utility::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Variantcalling::Allele_frequency qw{ calculate_af max_af };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view bcftools_index };
    use MIP::Program::Variantcalling::Vt qw{ vt_decompose vt_normalize vt_uniq };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP recipe name
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    my $file_path;
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
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $infile_path,
                output_type => $bcftools_output_type,
            }
        );
    }
    ## Replace #FORMAT field prior to smart decomposition (variant vcfs)
    if ($gnu_sed) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        gnu_sed(
            {
                FILEHANDLE => $FILEHANDLE,
                script     => q{'s/ID=AD,Number=./ID=AD,Number=R/'},
            }
        );
    }
    if ($decompose) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        vt_decompose(
            {
                FILEHANDLE             => $FILEHANDLE,
                infile_path            => q{-},
                smart_decomposition    => 1,
                stderrfile_path_append => $stderrfile_path,
            }
        );
    }
    if ($normalize) {

        # Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        vt_normalize(
            {
                FILEHANDLE                     => $FILEHANDLE,
                infile_path                    => q{-},
                no_fail_inconsistent_reference => 1,
                referencefile_path             => $human_genome_reference,
                stderrfile_path_append         => $stderrfile_path,
            }
        );
    }
    if ( $active_parameter_href->{vt_uniq} > 0 && $uniq ) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        vt_uniq(
            {
                FILEHANDLE             => $FILEHANDLE,
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
        print {$FILEHANDLE} $PIPE . $SPACE;

        htslib_bgzip(
            {
                FILEHANDLE      => $FILEHANDLE,
                write_to_stdout => 1,
            }
        );
    }

    ## Temporary outfile
    print {$FILEHANDLE} q{>}
      . $SPACE
      . $outfile_path
      . $UNDERSCORE
      . q{splitted}
      . $UNDERSCORE
      . $random_integer
      . $SPACE;
    print {$FILEHANDLE} $cmd_break;

    ## tabix index
    if ( -e $infile_path_tbi || $tabix ) {

        my $tabix_infile_path =
          $outfile_path . $UNDERSCORE . q{splitted} . $UNDERSCORE . $random_integer;
        htslib_tabix(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => $tabix_infile_path,
                preset      => q{vcf},
            }
        );
        print {$FILEHANDLE} $cmd_break;

        ## Move index in place
        gnu_mv(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $tabix_infile_path . $tabix_suffix,
                outfile_path => $outfile_path . $tabix_suffix,
            }
        );
        print {$FILEHANDLE} $cmd_break;
    }

    ## Move processed reference to original place
    gnu_mv(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $outfile_path
              . $UNDERSCORE
              . q{splitted}
              . $UNDERSCORE
              . $random_integer,
            outfile_path => $outfile_path,
        }
    );
    print {$FILEHANDLE} $cmd_break;

    return;
}

1;
