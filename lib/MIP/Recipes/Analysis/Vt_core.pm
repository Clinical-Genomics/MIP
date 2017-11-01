package MIP::Recipes::Analysis::Vt_core;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

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

## analysis_vt_core

## Function : Split multi allelic records into single records and normalize
## Returns  : ""
## Arguments: $parameter_href, $active_parameter_href, $infile_lane_prefix_href, $job_id_href, $infile_path, $FILEHANDLE, $contig, $family_id, $human_genome_reference, $outfile_path, $bcftools_output_type, $core_number, $decompose, $normalize, $uniq, $gnu_sed, $program_name, $program_directory, $bgzip, $tabix, $instream, $cmd_break, $xargs_file_path_prefix
##          : $parameter_href          => Hash with paremters from yaml file {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $infile_path             => Infile path
##          : $FILEHANDLE              => Filehandle to write to
##          : $contig                  => The contig to extract {OPTIONAL, REF}
##          : $family_id               => The family ID
##          : $human_genome_reference  => Human genome reference
##          : $outfile_path            => Outfile path
##          : $bcftools_output_type    => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $core_number             => The number of cores to allocate
##          : $decompose               => Vt program decomnpose for splitting multiallelic variants
##          : $normalize               => Vt program normalize for normalizing to reference used in analysis
##          : $uniq                    => Vt program uniq for removing variant duplication that appear later in file
##          : $gnu_sed                 => Sed program for changing vcf #FORMAT field in variant vcfs
##          : $program_name            => Program name
##          : $program_directory       => Program directory to write to in sbatch script
##          : $bgzip                   => Compress output from vt using bgzip
##          : $tabix                   => Index compressed output using tabix
##          : $instream                => Data to vt is supplied as a unix pipe
##          : $cmd_break               => Command line separator ['"\n\n"'|";"]
##          : $xargs_file_path_prefix  => The xargs sbatch script file name {OPTIONAL}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $human_genome_reference;
    my $outfile_path;
    my $bcftools_output_type;
    my $core_number;
    my $decompose;
    my $normalize;
    my $uniq;
    my $gnu_sed;
    my $program_name;
    my $program_directory;
    my $bgzip;
    my $tabix;
    my $instream;
    my $cmd_break;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $infile_path;
    my $FILEHANDLE;
    my $contig;
    my $xargs_file_path_prefix;

    my $tmpl = {
        parameter_href =>
          { default => {}, strict_type => 1, store => \$parameter_href },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        job_id_href =>
          { default => {}, strict_type => 1, store => \$job_id_href },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        FILEHANDLE => { store => \$FILEHANDLE },
        xargs_file_path_prefix =>
          { strict_type => 1, store => \$xargs_file_path_prefix },
        contig    => { strict_type => 1, store => \$contig },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        human_genome_reference => {
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference
        },
        ## Use same path as infile path unless parameter is supplied
        outfile_path => {
            default     => $arg_href->{infile_path},
            strict_type => 1,
            store       => \$outfile_path,
        },
        bcftools_output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$bcftools_output_type,
        },
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
        normalize => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$normalize
        },
        uniq => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$uniq
        },
        gnu_sed => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$gnu_sed
        },
        program_name =>
          { default => q{vt}, strict_type => 1, store => \$program_name },
        program_directory =>
          { default => q{vt}, strict_type => 1, store => \$program_directory },
        bgzip => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$bgzip
        },
        tabix => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$tabix
        },
        instream => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$instream
        },
        cmd_break =>
          { default => $NEWLINE x 2, strict_type => 1, store => \$cmd_break },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Gnu::Software::Gnu_sed qw{ gnu_sed };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_view bcftools_index };
    use Program::Variantcalling::Mip qw{ calculate_af max_af };
    use Program::Htslib qw{ bgzip tabix };
    use Program::Variantcalling::Vt qw{ decompose normalize vt_uniq };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_add_to_samples };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};

    my $file_path;
    my $program_info_path;

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ## Create anonymous filehandle
    $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $family_id,
            program_name          => $program_name,
            program_directory     => $program_directory,
            core_number           => $core_number,
            process_time          => 20,
        }
    );

    ### Split multi allelic records into single records and normalize

    ## Get parameters
    my $stderrfile_path;
    my $append_stderr_info;

    ## Write stderr for xargs process
    if ( defined $xargs_file_path_prefix && defined $contig ) {

        ## Redirect xargs output to program specific stderr file
        $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        $append_stderr_info = 1;
    }

    ## Initate processing
    if ( not $instream ) {

        bcftools_view(
            {
                infile_path => $infile_path,
                output_type => $bcftools_output_type,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
    }
    ## Replace #FORMAT field prior to smart decomposition (variant vcfs)
    if ($gnu_sed) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        gnu_sed(
            {
                script     => q{'s/ID=AD,Number=./ID=AD,Number=R/'},
                FILEHANDLE => $FILEHANDLE,
            }
        );
    }
    if ($decompose) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        decompose(
            {
                infile_path         => q{-},
                stderrfile_path     => $stderrfile_path,
                FILEHANDLE          => $FILEHANDLE,
                smart_decomposition => 1,
            }
        );
    }
    if ($normalize) {

        # Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        normalize(
            {
                infile_path                    => q{-},
                stderrfile_path                => $stderrfile_path,
                append_stderr_info             => 1,
                referencefile_path             => $human_genome_reference,
                no_fail_inconsistent_reference => 1,
                FILEHANDLE                     => $FILEHANDLE,
            }
        );
    }
    if ($uniq) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        vt_uniq(
            {
                infile_path        => q{-},
                stderrfile_path    => $stderrfile_path,
                append_stderr_info => $append_stderr_info,
                FILEHANDLE         => $FILEHANDLE,
            }
        );
    }

    ## tabix has been/will be used on file, compress again
    my $tabix_suffix    = $DOT . q{tbi};
    my $infile_path_tbi = $infile_path . $tabix_suffix;
    if ( -e $infile_path_tbi || $bgzip ) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        bgzip(
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
            $outfile_path
          . $UNDERSCORE
          . q{splitted}
          . $UNDERSCORE
          . $random_integer;
        tabix(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $tabix_infile_path,
                force       => 1,
                preset      => q{vcf},
            }
        );
        print {$FILEHANDLE} $cmd_break;

        ## Move index in place
        gnu_mv(
            {
                infile_path  => $tabix_infile_path . $tabix_suffix,
                outfile_path => $outfile_path . $tabix_suffix,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        print {$FILEHANDLE} $cmd_break;
    }

    ## Move processed reference to original place
    gnu_mv(
        {
            infile_path => $outfile_path
              . $UNDERSCORE
              . q{splitted}
              . $UNDERSCORE
              . $random_integer,
            outfile_path => $outfile_path,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    print {$FILEHANDLE} $cmd_break;

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

        slurm_submit_job_no_dependency_add_to_samples(
            {
                job_id_href      => $job_id_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $family_id,
                path             => $job_id_chain,
                sbatch_file_name => $file_path,
                log              => $log,
            }
        );
    }
    return;
}

sub analysis_vt_core_rio {

## analysis_vt_core

## Function : Split multi allelic records into single records and normalize
## Returns  : ""
## Arguments: $active_parameter_href, $infile_path, $FILEHANDLE, $contig, $family_id, $human_genome_reference, $outfile_path, $bcftools_output_type, $core_number, $decompose, $normalize, $uniq, $gnu_sed, $program_name, $program_directory, $bgzip, $tabix, $instream, $cmd_break, $xargs_file_path_prefix
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $infile_path             => Infile path
##          : $FILEHANDLE              => Filehandle to write to
##          : $contig                  => The contig to extract {OPTIONAL, REF}
##          : $family_id               => The family ID
##          : $human_genome_reference  => Human genome reference
##          : $outfile_path            => Outfile path
##          : $bcftools_output_type    => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $core_number             => The number of cores to allocate
##          : $decompose               => Vt program decomnpose for splitting multiallelic variants
##          : $normalize               => Vt program normalize for normalizing to reference used in analysis
##          : $uniq                    => Vt program uniq for removing variant duplication that appear later in file
##          : $gnu_sed                 => Sed program for changing vcf #FORMAT field in variant vcfs
##          : $program_name            => Program name
##          : $program_directory       => Program directory to write to in sbatch script
##          : $bgzip                   => Compress output from vt using bgzip
##          : $tabix                   => Index compressed output using tabix
##          : $instream                => Data to vt is supplied as a unix pipe
##          : $cmd_break               => Command line separator ['"\n\n"'|";"]
##          : $xargs_file_path_prefix  => The xargs sbatch script file name {OPTIONAL}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;
    my $human_genome_reference;
    my $outfile_path;
    my $bcftools_output_type;
    my $core_number;
    my $decompose;
    my $normalize;
    my $uniq;
    my $gnu_sed;
    my $program_name;
    my $program_directory;
    my $bgzip;
    my $tabix;
    my $instream;
    my $cmd_break;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $infile_path;
    my $FILEHANDLE;
    my $contig;
    my $xargs_file_path_prefix;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        FILEHANDLE => { store => \$FILEHANDLE },
        xargs_file_path_prefix =>
          { strict_type => 1, store => \$xargs_file_path_prefix },
        contig    => { strict_type => 1, store => \$contig },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id
        },
        human_genome_reference => {
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference
        },
        ## Use same path as infile path unless parameter is supplied
        outfile_path => {
            default     => $arg_href->{infile_path},
            strict_type => 1,
            store       => \$outfile_path,
        },
        bcftools_output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$bcftools_output_type,
        },
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
        normalize => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$normalize
        },
        uniq => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$uniq
        },
        gnu_sed => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$gnu_sed
        },
        program_name =>
          { default => q{vt}, strict_type => 1, store => \$program_name },
        program_directory =>
          { default => q{vt}, strict_type => 1, store => \$program_directory },
        bgzip => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$bgzip
        },
        tabix => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$tabix
        },
        instream => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$instream
        },
        cmd_break =>
          { default => $NEWLINE x 2, strict_type => 1, store => \$cmd_break },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Gnu::Software::Gnu_sed qw{ gnu_sed };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_view bcftools_index };
    use Program::Variantcalling::Mip qw{ calculate_af max_af };
    use Program::Htslib qw{ bgzip tabix };
    use Program::Variantcalling::Vt qw{ decompose normalize vt_uniq };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_add_to_samples };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    my $file_path;
    my $program_info_path;

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ### Split multi allelic records into single records and normalize

    ## Get parameters
    my $stderrfile_path;
    my $append_stderr_info;

    ## Write stderr for xargs process
    if ( defined $xargs_file_path_prefix && defined $contig ) {

        ## Redirect xargs output to program specific stderr file
        $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
        $append_stderr_info = 1;
    }

    ## Initate processing
    if ( not $instream ) {

        bcftools_view(
            {
                infile_path => $infile_path,
                output_type => $bcftools_output_type,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
    }
    ## Replace #FORMAT field prior to smart decomposition (variant vcfs)
    if ($gnu_sed) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        gnu_sed(
            {
                script     => q{'s/ID=AD,Number=./ID=AD,Number=R/'},
                FILEHANDLE => $FILEHANDLE,
            }
        );
    }
    if ($decompose) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        decompose(
            {
                infile_path         => q{-},
                stderrfile_path     => $stderrfile_path,
                FILEHANDLE          => $FILEHANDLE,
                smart_decomposition => 1,
            }
        );
    }
    if ($normalize) {

        # Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        normalize(
            {
                infile_path                    => q{-},
                stderrfile_path                => $stderrfile_path,
                append_stderr_info             => 1,
                referencefile_path             => $human_genome_reference,
                no_fail_inconsistent_reference => 1,
                FILEHANDLE                     => $FILEHANDLE,
            }
        );
    }
    if ( $active_parameter_href->{vt_uniq} > 0 && $uniq ) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        vt_uniq(
            {
                infile_path        => q{-},
                stderrfile_path    => $stderrfile_path,
                append_stderr_info => $append_stderr_info,
                FILEHANDLE         => $FILEHANDLE,
            }
        );
    }

    ## tabix has been/will be used on file, compress again
    my $tabix_suffix    = $DOT . q{tbi};
    my $infile_path_tbi = $infile_path . $tabix_suffix;
    if ( -e $infile_path_tbi || $bgzip ) {

        ## Pipe
        print {$FILEHANDLE} $PIPE . $SPACE;

        bgzip(
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
            $outfile_path
          . $UNDERSCORE
          . q{splitted}
          . $UNDERSCORE
          . $random_integer;
        tabix(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $tabix_infile_path,
                force       => 1,
                preset      => q{vcf},
            }
        );
        print {$FILEHANDLE} $cmd_break;

        ## Move index in place
        gnu_mv(
            {
                infile_path  => $tabix_infile_path . $tabix_suffix,
                outfile_path => $outfile_path . $tabix_suffix,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        print {$FILEHANDLE} $cmd_break;
    }

    ## Move processed reference to original place
    gnu_mv(
        {
            infile_path => $outfile_path
              . $UNDERSCORE
              . q{splitted}
              . $UNDERSCORE
              . $random_integer,
            outfile_path => $outfile_path,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    print {$FILEHANDLE} $cmd_break;

    return;
}

1;
