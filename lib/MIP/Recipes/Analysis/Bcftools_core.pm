package MIP::Recipes::Analysis::Bcftools_core;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
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

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_bcftools_core };
}

sub analysis_bcftools_core {

## Function : Split multi allelic records into single records and normalize
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $bcftools_output_type    => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $bgzip                   => Compress output from bcftools using bgzip
##          : $build_gatk_index        => Build gatk index (.idx) file
##          : $cmd_break               => Command line separator ['"\n\n"'|";"]
##          : $contig                  => The contig to extract {OPTIONAL, REF}
##          : $core_number             => The number of cores to allocate
##          : $case_id                 => The case ID
##          : $filehandle              => Filehandle to write to
##          : $human_genome_reference  => Human genome reference
##          : $job_id_href             => Job id hash {REF}
##          : $infile_path             => Infile path
##          : $outfile_path            => Outfile path
##          : $parameter_href          => Hash with paremters from yaml file {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_directory        => Program directory to write to in sbatch script
##          : $recipe_name             => Program name
##          : $tabix                   => Index compressed output using tabix
##          : $xargs_file_path_prefix  => The xargs sbatch script file name {OPTIONAL}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $contig;
    my $filehandle;
    my $job_id_href;
    my $infile_path;
    my $parameter_href;
    my $xargs_file_path_prefix;

    ## Default(s)
    my $bcftools_output_type;
    my $bgzip;
    my $build_gatk_index;
    my $cmd_break;
    my $core_number;
    my $case_id;
    my $human_genome_reference;
    my $outfile_path;
    my $profile_base_command;
    my $recipe_directory;
    my $recipe_name;
    my $tabix;

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
        build_gatk_index => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$build_gatk_index,
            strict_type => 1,
        },
        cmd_break   => { default => $NEWLINE x 2, store => \$cmd_break, strict_type => 1, },
        core_number => {
            allow       => qr/ \A \d+ \z /xsm,
            default     => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        contig  => { store => \$contig, strict_type => 1, },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        filehandle             => { store => \$filehandle, },
        human_genome_reference => {
            default     => $arg_href->{active_parameter_href}{human_genome_reference},
            store       => \$human_genome_reference,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        job_id_href  => { default => {}, store => \$job_id_href, strict_type => 1, },
        outfile_path => {
            default     => $arg_href->{infile_path},
            store       => \$outfile_path,
            strict_type => 1,
        },
        parameter_href       => { default => {}, strict_type => 1, store => \$parameter_href, },
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_directory =>
          { default => q{bcftools_core}, store => \$recipe_directory, strict_type => 1, },
        recipe_name => { default => q{bcftools_core}, store => \$recipe_name, strict_type => 1, },
        tabix       => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$tabix,
            strict_type => 1,
        },
        xargs_file_path_prefix => { store => \$xargs_file_path_prefix, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Coreutils qw{ gnu_mv };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Gatk qw{ gatk_indexfeaturefile };
    use MIP::Program::Htslib qw{ htslib_bgzip htslib_tabix };
    use MIP::Program::Bcftools qw{ bcftools_index bcftools_norm bcftools_view };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Set MIP recipe name
    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
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
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
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
        $stderrfile_path = $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};
    }

    ## Initate processing
    bcftools_view(
        {
            filehandle  => $filehandle,
            infile_path => $infile_path,
            output_type => $bcftools_output_type,
        }
    );

    ## Pipe
    print {$filehandle} $PIPE . $SPACE;

    bcftools_norm(
        {
            filehandle      => $filehandle,
            infile_path     => q{-},
            multiallelic    => q{-},
            reference_check => q{s},
            reference_path  => $human_genome_reference,
        }
    );

    ## Pipe
    print {$filehandle} $PIPE . $SPACE;

    bcftools_norm(
        {

            filehandle        => $filehandle,
            infile_path       => q{-},
            remove_duplicates => 1,
        }
    );

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

    if ($build_gatk_index) {

        gatk_indexfeaturefile(
            {
                filehandle           => $filehandle,
                infile_path          => $outfile_path,
                java_use_large_pages => $active_parameter_href->{java_use_large_pages},
                temp_directory       => $active_parameter_href->{temp_directory},
            }
        );
        say {$filehandle} $NEWLINE;
    }

    close $filehandle;

    if ( $recipe{mode} == 1 ) {

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

1;
