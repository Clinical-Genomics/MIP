package MIP::Recipes::Analysis::Peddy;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_peddy };

}

sub analysis_peddy {

## Function : Compares case-relationships and sexes.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Pedigree qw{ create_fam_file };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bcftools qw{ bcftools_view_and_index_vcf };
    use MIP::Program::Peddy qw{ peddy };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

    my $genome_reference_version = $file_info_href->{human_genome_reference_version};
    my %recipe                   = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my %peddy_outfile = (
        peddy     => q{ped},
        ped_check => q{csv},
        sex_check => q{csv},
    );
    my @peddy_outfiles =
      map { $_ . $DOT . $peddy_outfile{$_} } keys %peddy_outfile;
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@peddy_outfiles,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

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
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $recipe_info_file ) = splitpath($recipe_info_path);

    ### SHELL:

    my $case_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            case_id          => $case_id,
            fam_file_path    => $case_file_path,
            filehandle       => $filehandle,
            sample_ids_ref   => $active_parameter_href->{sample_ids},
            sample_info_href => $sample_info_href,
        }
    );

    ## Peddy
    my $genome_site = _get_genome_site( { version => $genome_reference_version } );
    peddy(
        {
            case_file_path      => $case_file_path,
            filehandle          => $filehandle,
            genome_site         => $genome_site,
            infile_path         => $infile_path,
            outfile_prefix_path => $outfile_path_prefix,
        }
    );
    say {$filehandle} $NEWLINE;

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

      PEDDY_OUTPUT_FILES:
        while ( my ( $outfile_tag, $outfile_path ) = each %outfile_path ) {

            ## Collect QC metadata info for later use
            set_recipe_metafile_in_sample_info(
                {
                    metafile_tag     => $outfile_tag,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );

            if ( $outfile_tag eq q{ped_check} ) {

                ## Duplicate ped_check tag one level out in sample_info. To be used for automatic kinship test
                set_recipe_outfile_in_sample_info(
                    {
                        path             => $outfile_path,
                        recipe_name      => $outfile_tag,
                        sample_info_href => $sample_info_href,
                    }
                );
            }

            set_file_path_to_store(
                {
                    format           => q{meta},
                    id               => $case_id,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                    tag              => $outfile_tag,
                }
            );
        }

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{case_to_island},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _get_genome_site {

## Function : Get the genomic site for peddy if using grch38
## Returns  : hg38
## Arguments: $version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $version;

    my $tmpl = {
        version => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 1,
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not $version eq q{38} );

    return q{hg38};
}
1;
