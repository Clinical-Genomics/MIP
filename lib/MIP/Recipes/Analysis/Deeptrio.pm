package MIP::Recipes::Analysis::Deeptrio;

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

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_deeptrio };

}

sub analysis_deeptrio {

## Function : Returns vcfs and gvcfs from bam files using DeepVariant's deeptrio caller
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
##          : $sample_info_href        => Info on samples and case hash {REF}

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Deeptrio qw{ deeptrio };
    use MIP::Sample_info
      qw{ get_case_members_attributes_in_duos get_family_member_id set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %case_members_id =
      scalar @{ $active_parameter_href->{sample_ids} } == 2
      ? get_case_members_attributes_in_duos( { sample_info_href => $sample_info_href } )
      : get_family_member_id( { sample_info_href => $sample_info_href } );

    ## Unpack parameters
    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );

    my $model_type = _consensus_analysis_type(
        {
            parameter_href => $parameter_href,
        }
    );

    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Get the io infiles per chain and id
    my $infile_name_prefix;
    my $output_gvcf_child;
    my $output_gvcf_parent1;
    my $output_gvcf_parent2;
    my $output_vcf_child;
    my $output_vcf_parent1;
    my $output_vcf_parent2;
    my $reads_child;
    my $reads_parent1;
    my $reads_parent2;
    my $sample_name_child;
    my $sample_name_parent1;
    my $sample_name_parent2;

    my @parents = _get_parents(
        {
            case_members_id_href => \%case_members_id,
        }
    );
    my $child_id = $case_members_id{children}[0];

  SAMPLE_ID:
    foreach my $sample_id ( @parents, $child_id ) {
        my %sample_bam_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => q{markduplicates},
                stream         => q{out},
            }
        );
        $infile_name_prefix = $sample_bam_io{in}{file_name_prefix};

        my %sample_vcf_io = (
            %sample_bam_io,
            parse_io_outfiles(
                {
                    chain_id               => $job_id_chain,
                    id                     => $sample_id,
                    file_info_href         => $file_info_href,
                    file_name_prefixes_ref => [$infile_name_prefix],
                    outdata_dir            => $active_parameter_href->{outdata_dir},
                    parameter_href         => $parameter_href,
                    recipe_name            => $recipe_name,
                }
            )
        );
        if ( not $sample_name_parent1 and $sample_id ne $child_id ) {
            $reads_parent1 =
              $sample_bam_io{out}{file_path_prefix} . $sample_bam_io{out}{file_suffix};
            $output_gvcf_parent1 = $sample_vcf_io{out}{file_path};
            $output_vcf_parent1  = $sample_vcf_io{out}{file_path_prefix} . q{.vcf.gz};
            $sample_name_parent1 = $sample_id;
        }
        elsif ( $sample_id ne $child_id ) {
            $reads_parent2 =
              $sample_bam_io{out}{file_path_prefix} . $sample_bam_io{out}{file_suffix};
            $output_gvcf_parent2 = $sample_vcf_io{out}{file_path};
            $output_vcf_parent2  = $sample_vcf_io{out}{file_path_prefix} . q{.vcf.gz};
            $sample_name_parent2 = $sample_id;
        }
        else {
            $reads_child = $sample_bam_io{out}{file_path_prefix} . $sample_bam_io{out}{file_suffix};
            $output_gvcf_child = $sample_vcf_io{out}{file_path};
            $output_vcf_child  = $sample_vcf_io{out}{file_path_prefix} . q{.vcf.gz};
            $sample_name_child = $sample_id;
        }
    }

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe_resource{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            gpu_number            => $recipe_resource{gpu_number},
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe_resource{memory},
            process_time          => $recipe_resource{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    deeptrio(
        {
            filehandle          => $filehandle,
            model_type          => $model_type,
            num_shards          => $recipe_resource{core_number},
            output_gvcf_child   => $output_gvcf_child,
            output_gvcf_parent1 => $output_gvcf_parent1,
            output_gvcf_parent2 => $output_gvcf_parent2,
            output_vcf_child    => $output_vcf_child,
            output_vcf_parent1  => $output_vcf_parent1,
            output_vcf_parent2  => $output_vcf_parent2,
            referencefile_path  => $active_parameter_href->{human_genome_reference},
            reads_child         => $reads_child,
            reads_parent1       => $reads_parent1,
            reads_parent2       => $reads_parent2,
            sample_name_child   => $sample_name_child,
            sample_name_parent1 => $sample_name_parent1,
            sample_name_parent2 => $sample_name_parent2,
        }
    );

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $output_gvcf_child,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $job_id_chain,
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

sub _consensus_analysis_type {

## Function : Determine model type from consensus analysis type
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %consensus_analysis_type_map =
      ( MIXED => q{WGS}, PANEL => q{WES}, WGS => q{WGS}, WES => q{WES} );

    my $consensus_analysis_type =
      $consensus_analysis_type_map{ uc $parameter_href->{cache}{consensus_analysis_type} };

    return $consensus_analysis_type;
}

sub _get_parents {

## Function : Determine model type from consensus analysis type
## Returns  :
## Arguments: $case_members_id => Case member id hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_members_id_href;

    my $tmpl = {
        case_members_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$case_members_id_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};
    my @parents;

    foreach my $parent (qw{ father mother }) {

        if ( exists $case_members_id_href->{$parent} ) {
            push @parents, $case_members_id_href->{$parent};
        }
    }
    return @parents;
}

1;
