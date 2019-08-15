package MIP::Main::Analyse;

#### Master script for analysing paired end reads from the Illumina plattform in fastq(.gz) format to annotated ranked disease causing variants. The program performs QC, aligns reads using BWA, performs variant discovery and annotation as well as ranking the found variants according to disease potential.

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ basename fileparse };
use File::Copy qw{ copy };
use File::Spec::Functions qw{ catdir catfile devnull };
use FindBin qw{ $Bin };
use Getopt::Long;
use IPC::Cmd qw{ can_run run};
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use POSIX;
use Time::Piece;
use utf8;
use warnings qw{ FATAL utf8 };

## Third party module(s)
use autodie qw{ open close :all };
use IPC::System::Simple;
use List::MoreUtils qw { any uniq all };
use Modern::Perl qw{ 2014 };
use Path::Iterator::Rule;
use Readonly;

## MIPs lib/
# Add MIPs internal lib
use MIP::Check::Modules qw{ check_perl_modules };
use MIP::Check::Parameter qw{ check_allowed_temp_directory
  check_cmd_config_vs_definition_file
  check_email_address
  check_load_env_packages
  check_parameter_hash
  check_recipe_exists_in_hash
  check_recipe_name
  check_recipe_mode
  check_sample_ids
};
use MIP::Check::Path qw{ check_executable_in_path check_parameter_files };
use MIP::Check::Reference qw{ check_human_genome_file_endings };
use MIP::Constants qw{ $DOT $EMPTY_STR $MIP_VERSION $NEWLINE $SINGLE_QUOTE $SPACE $TAB };
use MIP::Cluster qw{ check_max_core_number check_recipe_memory_allocation };
use MIP::File::Format::Mip qw{ build_file_prefix_tag };
use MIP::File::Format::Pedigree
  qw{ create_fam_file detect_founders detect_sample_id_gender detect_trio parse_yaml_pedigree_file reload_previous_pedigree_info };
use MIP::File::Format::Yaml qw{ load_yaml write_yaml order_parameter_names };
use MIP::Get::Analysis qw{ get_overall_analysis_type };
use MIP::Get::Parameter qw{ get_program_executables };
use MIP::Log::MIP_log4perl qw{ initiate_logger set_default_log4perl_file };
use MIP::Parse::Parameter qw{ parse_start_with_recipe };
use MIP::Processmanagement::Processes qw{ write_job_ids_to_file };
use MIP::Set::Contigs qw{ set_contigs };
use MIP::Set::Parameter
  qw{ set_config_to_active_parameters set_custom_default_to_active_parameter set_default_config_dynamic_parameters set_default_to_active_parameter set_cache set_human_genome_reference_features set_no_dry_run_parameters set_parameter_reference_dir_path set_recipe_resource };
use MIP::Update::Parameters
  qw{ update_dynamic_config_parameters update_reference_parameters update_vcfparser_outfile_counter };
use MIP::Update::Path qw{ update_to_absolute_path };
use MIP::Update::Recipes qw{ update_recipe_mode_with_dry_run_all };

## Recipes
use MIP::Recipes::Pipeline::Analyse_dragen_rd_dna qw{ pipeline_analyse_dragen_rd_dna };
use MIP::Recipes::Pipeline::Analyse_rd_dna qw{ pipeline_analyse_rd_dna };
use MIP::Recipes::Pipeline::Analyse_rd_rna qw{ pipeline_analyse_rd_rna };
use MIP::Recipes::Pipeline::Analyse_rd_dna_vcf_rerun
  qw{ pipeline_analyse_rd_dna_vcf_rerun };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.17;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_analyse };
}

sub mip_analyse {

## Function : Execute mip analyse pre pipeline parsing
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href        => File info hash {REF}
##          : $order_parameters_ref  => Order of addition to parameter array {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $order_parameters_ref;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Transfer to lexical variables
    # Parameters to include in each analysis run
    my %active_parameter = %{$active_parameter_href};

    # Holds all active parameters values for broadcasting
    my @broadcasts;

    # File information
    my %file_info = %{$file_info_href};

    # Order parameters for logical broadcast of parameters
    my @order_parameters = @{$order_parameters_ref};

    # All parameters MIP analyse knows
    my %parameter = %{$parameter_href};

#### Script parameters

## Add date_time_stamp for later use in log and qc_metrics yaml file
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;
    my $date            = $date_time->ymd;

    # Catches script name and removes ending
    my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{pl} ) );
    chomp( $date_time_stamp, $date, $script );

#### Set program parameters

## Directories, files, job_ids and sample_info
    my ( %infile_lane_prefix, %infile_both_strands_prefix, %job_id, %sample_info );

#### Staging Area
### Get and/or set input parameters

## Change relative path to absolute path for parameter with "update_path: absolute_path" in config
    update_to_absolute_path(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );

### Config file
## If config from cmd
    if ( exists $active_parameter{config_file}
        && defined $active_parameter{config_file} )
    {

        ## Loads a YAML file into an arbitrary hash and returns it.
        my %config_parameter =
          load_yaml( { yaml_file => $active_parameter{config_file}, } );

        ## Remove previous analysis specific info not relevant for current run e.g. log file, which is read from pedigree or cmd
        my @remove_keys = (
            qw{ found_female found_male found_other found_other_count log_file dry_run_all }
        );

      KEY:
        foreach my $key (@remove_keys) {

            delete $config_parameter{$key};
        }

## Set config parameters into %active_parameter unless $parameter
## has been supplied on the command line
        set_config_to_active_parameters(
            {
                active_parameter_href => \%active_parameter,
                config_parameter_href => \%config_parameter,
            }
        );

        ## Compare keys from config and cmd (%active_parameter) with definitions file (%parameter)
        check_cmd_config_vs_definition_file(
            {
                active_parameter_href => \%active_parameter,
                parameter_href        => \%parameter,
            }
        );

        my @config_dynamic_parameters = qw{ analysis_constant_path };

        ## Replace config parameter with cmd info for config dynamic parameter
        set_default_config_dynamic_parameters(
            {
                active_parameter_href => \%active_parameter,
                parameter_href        => \%parameter,
                parameter_names_ref   => \@config_dynamic_parameters,
            }
        );

        ## Loop through all parameters and update info
      PARAMETER:
        foreach my $parameter_name ( keys %parameter ) {

            ## Updates the active parameters to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
            update_dynamic_config_parameters(
                {
                    active_parameter_href => \%active_parameter,
                    parameter_name        => $parameter_name,
                }
            );
        }
    }

## Set the default Log4perl file using supplied dynamic parameters.
    $active_parameter{log_file} = set_default_log4perl_file(
        {
            cmd_input       => $active_parameter{log_file},
            date            => $date,
            date_time_stamp => $date_time_stamp,
            outdata_dir     => $active_parameter{outdata_dir},
            script          => $script,
        }
    );

## Creates log object
    my $log = initiate_logger(
        {
            file_path => $active_parameter{log_file},
            log_name  => uc q{mip_analyse},
        }
    );

## Write MIP VERSION and log file path
    $log->info( q{MIP Version: } . $MIP_VERSION );
    $log->info( q{Script parameters and info from are saved in file: }
          . $active_parameter{log_file} );

## Parse pedigree file
## Reads case_id_pedigree file in YAML format. Checks for pedigree data for allowed entries and correct format. Add data to sample_info depending on user info.
    # Meta data in YAML format
    if ( defined $active_parameter{pedigree_file} ) {

        ## Loads a YAML file into an arbitrary hash and returns it. Load parameters from previous run from sample_info_file
        my %pedigree =
          load_yaml( { yaml_file => $active_parameter{pedigree_file}, } );

        $log->info( q{Loaded: } . $active_parameter{pedigree_file} );

        parse_yaml_pedigree_file(
            {
                active_parameter_href => \%active_parameter,
                file_path             => $active_parameter{pedigree_file},
                log                   => $log,
                parameter_href        => \%parameter,
                pedigree_href         => \%pedigree,
                sample_info_href      => \%sample_info,
            }
        );
    }

    # Detect if all samples has the same sequencing type and return consensus if reached
    $parameter{cache}{consensus_analysis_type} = get_overall_analysis_type(
        {
            analysis_type_href => \%{ $active_parameter{analysis_type} },
            log                => $log,
        }
    );

### Populate uninitilized active_parameters{parameter_name} with default from parameter
  PARAMETER:
    foreach my $parameter_name ( keys %parameter ) {

        ## If hash and set - skip
        next PARAMETER
          if ( ref $active_parameter{$parameter_name} eq qw{HASH}
            && keys %{ $active_parameter{$parameter_name} } );

        ## If array and set - skip
        next PARAMETER
          if ( ref $active_parameter{$parameter_name} eq qw{ARRAY}
            && @{ $active_parameter{$parameter_name} } );

        ## If scalar and set - skip
        next PARAMETER
          if ( defined $active_parameter{$parameter_name}
            and not ref $active_parameter{$parameter_name} );

        ### Special case for parameters that are dependent on other parameters values
        my @custom_default_parameters = qw{
          analysis_type
          bwa_build_reference
          exome_target_bed
          fusion_filter_reference_genome
          gatk_path
          infile_dirs
          picardtools_path
          salmon_quant_reference_genome
          sample_info_file
          snpeff_path
          star_aln_reference_genome
          rtg_vcfeval_reference_genome
          temp_directory
          vep_directory_path
        };

        if ( any { $_ eq $parameter_name } @custom_default_parameters ) {

            set_custom_default_to_active_parameter(
                {
                    active_parameter_href => \%active_parameter,
                    parameter_href        => \%parameter,
                    parameter_name        => $parameter_name,
                }
            );
            next PARAMETER;
        }

        ## Checks and sets user input or default values to active_parameters
        set_default_to_active_parameter(
            {
                active_parameter_href => \%active_parameter,
                associated_recipes_ref =>
                  \@{ $parameter{$parameter_name}{associated_recipe} },
                parameter_href => \%parameter,
                parameter_name => $parameter_name,
            }
        );
    }

## Update path for supplied reference(s) associated with parameter that should reside in the mip reference directory to full path
    set_parameter_reference_dir_path(
        {
            active_parameter_href => \%active_parameter,
            parameter_name        => q{human_genome_reference},
        }
    );

## Detect version and source of the human_genome_reference: Source (hg19 or GRCh) and check compression status
    set_human_genome_reference_features(
        {
            file_info_href => \%file_info,
            human_genome_reference =>
              basename( $active_parameter{human_genome_reference} ),
            log            => $log,
            parameter_href => \%parameter,
        }
    );

## Reference in MIP reference directory
  PARAMETER:
    foreach my $parameter_name ( keys %parameter ) {

        ## Expect file to be in reference directory
        if ( exists $parameter{$parameter_name}{reference} ) {

            update_reference_parameters(
                {
                    active_parameter_href => \%active_parameter,
                    associated_recipes_ref =>
                      \@{ $parameter{$parameter_name}{associated_recipe} },
                    parameter_name => $parameter_name,
                }
            );
        }
    }

### Checks

## Check existence of files and directories
  PARAMETER:
    foreach my $parameter_name ( keys %parameter ) {

        if ( exists $parameter{$parameter_name}{exists_check} ) {

            check_parameter_files(
                {
                    active_parameter_href => \%active_parameter,
                    associated_recipes_ref =>
                      \@{ $parameter{$parameter_name}{associated_recipe} },
                    log                    => $log,
                    parameter_exists_check => $parameter{$parameter_name}{exists_check},
                    parameter_href         => \%parameter,
                    parameter_name         => $parameter_name,
                }
            );
        }
    }

## Updates sample_info hash with previous run pedigree info
    reload_previous_pedigree_info(
        {
            log                   => $log,
            sample_info_href      => \%sample_info,
            sample_info_file_path => $active_parameter{sample_info_file},
        }
    );

## Special case since dict is created with .fastq removed
## Check the existance of associated human genome files
    check_human_genome_file_endings(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            log                   => $log,
            parameter_href        => \%parameter,
            parameter_name        => q{human_genome_reference_file_endings},
        }
    );

## Detect case constellation based on pedigree file
    $parameter{cache}{trio} = detect_trio(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            sample_info_href      => \%sample_info,
        }
    );

## Detect number of founders (i.e. parents ) based on pedigree file
    detect_founders(
        {
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
        }
    );

## Check email adress syntax and mail host
    check_email_address(
        {
            email => $active_parameter{email},
            log   => $log,
        }
    );

## Check that the temp directory value is allowed
    check_allowed_temp_directory(
        {
            log            => $log,
            temp_directory => $active_parameter{temp_directory},
        }
    );

## Parameters that have keys as MIP recipe names
    my @parameter_keys_to_check = (
        qw{ recipe_time recipe_core_number recipe_memory
          set_recipe_core_number set_recipe_memory set_recipe_time }
    );
  PARAMETER_NAME:
    foreach my $parameter_name (@parameter_keys_to_check) {

        ## Test if key from query hash exists truth hash
        check_recipe_exists_in_hash(
            {
                log            => $log,
                parameter_name => $parameter_name,
                query_ref      => \%{ $active_parameter{$parameter_name} },
                truth_href     => \%parameter,
            }
        );
    }

## Set recipe resource allocation for specific recipe(s)
    set_recipe_resource( { active_parameter_href => \%active_parameter, } );

## Parameters with key(s) that have elements as MIP recipe names
    my @parameter_element_to_check = qw{ associated_recipe };
  PARAMETER:
    foreach my $parameter ( keys %parameter ) {

      KEY:
        foreach my $parameter_name (@parameter_element_to_check) {

            next KEY if ( not exists $parameter{$parameter}{$parameter_name} );

            ## Test if element from query array exists truth hash
            check_recipe_exists_in_hash(
                {
                    log            => $log,
                    parameter_name => $parameter_name,
                    query_ref      => \@{ $parameter{$parameter}{$parameter_name} },
                    truth_href     => \%parameter,
                }
            );
        }
    }

## Parameters that have elements as MIP recipe names
    my @parameter_elements_to_check =
      (qw(associated_recipe decompose_normalize_references));
    foreach my $parameter_name (@parameter_elements_to_check) {

        ## Test if element from query array exists truth hash
        check_recipe_exists_in_hash(
            {
                log            => $log,
                parameter_name => $parameter_name,
                query_ref      => \@{ $active_parameter{$parameter_name} },
                truth_href     => \%parameter,
            }
        );
    }

## Check that the recipe core number do not exceed the maximum per node
    foreach my $recipe_name ( keys %{ $active_parameter{recipe_core_number} } ) {

        ## Limit number of cores requested to the maximum number of cores available per node
        $active_parameter{recipe_core_number}{$recipe_name} = check_max_core_number(
            {
                max_cores_per_node => $active_parameter{max_cores_per_node},
                core_number_requested =>
                  $active_parameter{recipe_core_number}{$recipe_name},
            }
        );
    }

    ## Check that the recipe memory do not exceed the maximum per node
    foreach my $recipe_name ( keys %{ $active_parameter{recipe_memory} } ) {

        check_recipe_memory_allocation(
            {
                node_ram_memory => $active_parameter{node_ram_memory},
                recipe_memory_allocation =>
                  $active_parameter{recipe_memory}{$recipe_name},
            }
        );
    }

    ## Check programs in path, and executable
    check_executable_in_path(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    );

## Test that the case_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
    check_sample_ids(
        {
            case_id        => $active_parameter{case_id},
            log            => $log,
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    );

## Adds dynamic aggregate information from definitions to parameter hash
    set_cache(
        {
            aggregates_ref => [
                ## Collect all aligners
                q{recipe_type:aligners},
                ## Collects all references in that are supposed to be in reference directory
                q{reference:reference_dir},
                ## Collects all structural variant_callers
                q{recipe_type:structural_variant_callers},
                ## Collects all variant_callers
                q{recipe_type:variant_callers},
                ## Collects all recipes that MIP can handle
                q{type:recipe},
            ],
            parameter_href => \%parameter,
        }
    );

    @{ $parameter{cache}{program_executables} } =
      get_program_executables( { parameter_href => \%parameter, } );

## Check correct value for recipe mode in MIP
    check_recipe_mode(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    );

    ## Check that package name name are included in MIP as either "mip", "recipe" or "program"
    check_load_env_packages(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );

## Check that recipe name and program name are not identical
    check_recipe_name(
        {
            parameter_href   => \%parameter,
            recipe_names_ref => \@{ $parameter{cache}{recipe} },
        }
    );

    parse_start_with_recipe(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        },
    );

## Update recipe mode depending on dry_run_all flag
    update_recipe_mode_with_dry_run_all(
        {
            active_parameter_href => \%active_parameter,
            dry_run_all           => $active_parameter{dry_run_all},
            recipes_ref           => \@{ $parameter{cache}{recipe} },
        }
    );

## Detect the gender(s) included in current analysis
    (

        $active_parameter{found_male},
        $active_parameter{found_female},
        $active_parameter{found_other},
        $active_parameter{found_other_count},
      )
      = detect_sample_id_gender(
        {
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
        }
      );

### Contigs
## Set contig prefix and contig names depending on reference used
    set_contigs(
        {
            file_info_href => \%file_info,
            version        => $file_info{human_genome_reference_version},
        }
    );

## Creates all fileendings as the samples is processed depending on the chain of modules activated
    build_file_prefix_tag(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            order_recipes_ref     => \@{ $parameter{cache}{order_recipes_ref} },
            parameter_href        => \%parameter,
        }
    );

## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => \%active_parameter,
            execution_mode        => q{system},
            fam_file_path         => catfile(
                $active_parameter{outdata_dir}, $active_parameter{case_id},
                $active_parameter{case_id} . $DOT . q{fam}
            ),
            log              => $log,
            parameter_href   => \%parameter,
            sample_info_href => \%sample_info,
        }
    );

############
####MAIN####
############

    set_no_dry_run_parameters(
        {
            is_dry_run_all   => $active_parameter{dry_run_all},
            analysis_date    => $date_time_stamp,
            mip_version      => $MIP_VERSION,
            sample_info_href => \%sample_info,
        }
    );

    my $consensus_analysis_type = $parameter{cache}{consensus_analysis_type};

    ## Create dispatch table of pipelines
    my %pipeline = (
        dragen_rd_dna => \&pipeline_analyse_dragen_rd_dna,
        mixed         => \&pipeline_analyse_rd_dna,
        vrn           => \&pipeline_analyse_rd_dna_vcf_rerun,
        wes           => \&pipeline_analyse_rd_dna,
        wgs           => \&pipeline_analyse_rd_dna,
        wts           => \&pipeline_analyse_rd_rna,
    );

    $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );
    $pipeline{$consensus_analysis_type}->(
        {
            active_parameter_href           => \%active_parameter,
            broadcasts_ref                  => \@broadcasts,
            file_info_href                  => \%file_info,
            infile_both_strands_prefix_href => \%infile_both_strands_prefix,
            infile_lane_prefix_href         => \%infile_lane_prefix,
            job_id_href                     => \%job_id,
            log                             => $log,
            order_parameters_ref            => \@order_parameters,
            order_recipes_ref               => \@{ $parameter{cache}{order_recipes_ref} },
            parameter_href                  => \%parameter,
            sample_info_href                => \%sample_info,
        }
    );

## Write QC for recipes used in analysis
    # Write sample info to yaml file
    if ( $active_parameter{sample_info_file} ) {

        ## Writes a YAML hash to file
        write_yaml(
            {
                yaml_href      => \%sample_info,
                yaml_file_path => $active_parameter{sample_info_file},
            }
        );
        $log->info( q{Wrote: } . $active_parameter{sample_info_file} );
    }

## Write job_ids to file
    write_job_ids_to_file(
        {
            active_parameter_href => \%active_parameter,
            date_time_stamp       => $date_time_stamp,
            job_id_href           => \%job_id,
        }
    );

    return;
}

######################
####Sub routines######
######################

##Investigate potential autodie error
if ( $EVAL_ERROR and $EVAL_ERROR->isa(q{autodie::exception}) ) {

    if ( $EVAL_ERROR->matches(q{default}) ) {

        say {*STDERR} q{Not an autodie error at all};
    }
    if ( $EVAL_ERROR->matches(q{open}) ) {

        say {*STDERR} q{Error from open};
    }
    if ( $EVAL_ERROR->matches(q{:io}) ) {

        say {*STDERR} q{Non-open, IO error.};
    }
}
elsif ($EVAL_ERROR) {

    say {*STDERR} q{A non-autodie exception.};
}

1;
