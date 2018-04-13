package MIP::Main::Analyse;

#### Master script for analysing paired end reads from the Illumina plattform in fastq(.gz) format to annotated ranked disease causing variants. The program performs QC, aligns reads using BWA, performs variant discovery and annotation as well as ranking the found variants according to disease potential.

#### Copyright 2011 Henrik Stranneheim

use 5.018;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename fileparse };
use File::Copy qw{ copy };
use File::Path qw{ make_path };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
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
use MIP::Check::Cluster qw{ check_max_core_number };
use MIP::Check::Modules qw{ check_perl_modules };
use MIP::Check::Parameter
  qw{ check_allowed_temp_directory check_cmd_config_vs_definition_file check_email_address check_parameter_hash };
use MIP::Check::Path qw{ check_target_bed_file_suffix check_parameter_files };
use MIP::Check::Reference
  qw{ check_human_genome_file_endings check_parameter_metafiles };
use MIP::File::Format::Pedigree
  qw{ create_fam_file parse_yaml_pedigree_file reload_previous_pedigree_info };
use MIP::File::Format::Yaml qw{ load_yaml write_yaml order_parameter_names };
use MIP::Get::Analysis qw{ get_dependency_tree get_overall_analysis_type };
use MIP::Get::File qw{ get_select_file_contigs };
use MIP::Log::MIP_log4perl qw{ initiate_logger set_default_log4perl_file };
use MIP::Script::Utils qw{ help };
use MIP::Set::Contigs qw{ set_contigs };
use MIP::Set::Parameter
  qw{ set_config_to_active_parameters set_custom_default_to_active_parameter set_default_config_dynamic_parameters set_default_to_active_parameter set_dynamic_parameter set_human_genome_reference_features set_parameter_reference_dir_path set_parameter_to_broadcast };
use MIP::Update::Contigs qw{ update_contigs_for_run };
use MIP::Update::Parameters
  qw{ update_dynamic_config_parameters update_exome_target_bed update_reference_parameters update_vcfparser_outfile_counter };
use MIP::Update::Path qw{ update_to_absolute_path };
use MIP::Update::Programs
  qw{ update_prioritize_flag update_program_mode_for_analysis_type update_program_mode_with_dry_run_all update_program_mode_with_start_with };

## Recipes
use MIP::Recipes::Analysis::Gzip_fastq qw{ analysis_gzip_fastq };
use MIP::Recipes::Analysis::Split_fastq_file qw{ analysis_split_fastq_file };
use MIP::Recipes::Analysis::Vt_core qw{ analysis_vt_core };
use MIP::Recipes::Pipeline::Rare_disease qw{ pipeline_rare_disease };
use MIP::Recipes::Pipeline::Rna qw{ pipeline_rna };
use MIP::Recipes::Pipeline::Cancer qw{ pipeline_cancer };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_analyse };
}

## Constants
Readonly my $DOT       => q{.};
Readonly my $EMPTY_STR => q{};
Readonly my $NEWLINE   => qq{\n};
Readonly my $SPACE     => q{ };
Readonly my $TAB       => qq{\t};

sub mip_analyse {

## Function : Creates program directories (info & programData & programScript), program script filenames and writes sbatch header.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href                       => File info hash {REF}
##          : $parameter_href        => Parameter hash {REF}
#           : $order_parameters_ref  => Order of addition to parameter array {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $file_info_href;
    my $order_parameters_ref;

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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Transfer to lexical variables
    my %active_parameter = %{$active_parameter_href};
    my %file_info        = %{$file_info_href};
    my @order_parameters = @{$order_parameters_ref};
    my %parameter        = %{$parameter_href};

#### Script parameters

## Add date_time_stamp for later use in log and qc_metrics yaml file
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;
    my $date            = $date_time->ymd;

    # Catches script name and removes ending
    my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{pl} ) );
    chomp( $date_time_stamp, $date, $script );

#### Set program parameters

## Set MIP version
    our $VERSION = 'v7.0.1';

    if ( $active_parameter{version} ) {

        say STDOUT $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION,
          $NEWLINE;
        exit;
    }

## Directories, files, job_ids and sample_info
    my ( %infile, %indir_path, %infile_lane_prefix, %lane,
        %infile_both_strands_prefix, %job_id, %sample_info );

#### Staging Area
### Get and/or set input parameters

## Special case for boolean flag that will be removed from
## config upon loading
    my @boolean_parameter = qw{dry_run_all};
    foreach my $parameter (@boolean_parameter) {

        if ( not defined $active_parameter{$parameter} ) {

            delete $active_parameter{$parameter};
        }
    }

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
        my @remove_keys = (qw{ log_file dry_run_all });

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

        my @config_dynamic_parameters =
          qw{ analysis_constant_path outaligner_dir };

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
        foreach my $parameter_name (@order_parameters) {

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
            active_parameter_href => \%active_parameter,
            cmd_input             => $active_parameter{log_file},
            date                  => $date,
            date_time_stamp       => $date_time_stamp,
            script                => $script,
        }
    );

## Creates log object
    my $log = initiate_logger(
        {
            file_path => $active_parameter{log_file},
            log_name  => q{MIP},
        }
    );

## Parse pedigree file
## Reads family_id_pedigree file in YAML format. Checks for pedigree data for allowed entries and correct format. Add data to sample_info depending on user info.
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
                parameter_href        => \%parameter,
                pedigree_href         => \%pedigree,
                sample_info_href      => \%sample_info,
            }
        );
    }

# Detect if all samples has the same sequencing type and return consensus if reached
    $parameter{dynamic_parameter}{consensus_analysis_type} =
      get_overall_analysis_type(
        { analysis_type_href => \%{ $active_parameter{analysis_type} }, } );

### Populate uninitilized active_parameters{parameter_name} with default from parameter
  PARAMETER:
    foreach my $parameter_name (@order_parameters) {

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
        my @custom_default_parameters =
          qw{ analysis_type bwa_build_reference exome_target_bed infile_dirs sample_info_file rtg_vcfeval_reference_genome };

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
                associated_programs_ref =>
                  \@{ $parameter{$parameter_name}{associated_program} },
                log            => $log,
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

## Detect version and source of the human_genome_reference: Source (hg19 or GRCh).
    set_human_genome_reference_features(
        {
            file_info_href => \%file_info,
            human_genome_reference =>
              basename( $active_parameter{human_genome_reference} ),
            log => $log,
        }
    );

## Update exome_target_bed files with human_genome_reference_source and human_genome_reference_version
    update_exome_target_bed(
        {
            exome_target_bed_file_href => $active_parameter{exome_target_bed},
            human_genome_reference_source =>
              $file_info{human_genome_reference_source},
            human_genome_reference_version =>
              $file_info{human_genome_reference_version},
        }
    );

    # Holds all active parameters values for broadcasting
    my @broadcasts;

    if ( $active_parameter{verbose} ) {

        set_parameter_to_broadcast(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                order_parameters_ref  => \@order_parameters,
                broadcasts_ref        => \@broadcasts,
            }
        );
    }

## Reference in MIP reference directory
  PARAMETER:
    foreach my $parameter_name ( keys %parameter ) {

        ## Expect file to be in reference directory
        if ( exists $parameter{$parameter_name}{reference} ) {

            update_reference_parameters(
                {
                    active_parameter_href => \%active_parameter,
                    associated_programs_ref =>
                      \@{ $parameter{$parameter_name}{associated_program} },
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
                    associated_programs_ref =>
                      \@{ $parameter{$parameter_name}{associated_program} },
                    log => $log,
                    parameter_exists_check =>
                      $parameter{$parameter_name}{exists_check},
                    parameter_href => \%parameter,
                    parameter_name => $parameter_name,
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
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            parameter_name        => q{human_genome_reference},
        }
    );

## Check that supplied target file ends with ".bed" and otherwise croaks
  TARGET_FILE:
    foreach
      my $target_bed_file ( keys %{ $active_parameter{exome_target_bed} } )
    {

        check_target_bed_file_suffix(
            {
                parameter_name => q{exome_target_bed},
                path           => $target_bed_file,
            }
        );
    }

## Checks parameter metafile exists and set build_file parameter
    check_parameter_metafiles(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
        }
    );

## Update the expected number of outfile after vcfparser
    update_vcfparser_outfile_counter(
        { active_parameter_href => \%active_parameter, } );

## Collect select file contigs to loop over downstream
    if ( $active_parameter{vcfparser_select_file} ) {

## Collects sequences contigs used in select file
        @{ $file_info{select_file_contigs} } = get_select_file_contigs(
            {
                select_file_path =>
                  catfile( $active_parameter{vcfparser_select_file} ),
                log => $log,
            }
        );
    }

## Detect family constellation based on pedigree file
    $parameter{dynamic_parameter}{trio} = detect_trio(
        {
            active_parameter_href => \%active_parameter,
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
    if ( defined $active_parameter{email} ) {

        check_email_address(
            {
                email => $active_parameter{email},
                log   => $log,
            }
        );
    }

    if (
        not check_allowed_temp_directory(
            { temp_directory => $active_parameter{temp_directory}, }
        )
      )
    {

        $log->fatal( q{'--temp_directory }
              . $active_parameter{temp_directory}
              . q{' is not allowed because MIP will remove the temp directory after processing.}
              . "\n" );
        exit 1;
    }

## Parameters that have keys as MIP program names
    my @parameter_keys_to_check =
      (qw{ module_time module_core_number module_source_environment_command });
    foreach my $parameter_name (@parameter_keys_to_check) {

        ## Test if key from query hash exists truth hash
        check_key_exists_in_hash(
            {
                truth_href     => \%parameter,
                query_href     => \%{ $active_parameter{$parameter_name} },
                parameter_name => $parameter_name,
            }
        );
    }

## Check that the module core number do not exceed the maximum per node
    foreach my $program_name ( keys %{ $active_parameter{module_core_number} } )
    {

        ## Limit number of cores requested to the maximum number of cores available per node
        $active_parameter{module_core_number}{$program_name} =
          check_max_core_number(
            {
                max_cores_per_node => $active_parameter{max_cores_per_node},
                core_number_requested =>
                  $active_parameter{module_core_number}{$program_name},
            }
          );
    }

## Parameters that have elements as MIP program names
    my @parameter_elements_to_check =
      (qw(associated_program decompose_normalize_references));
    foreach my $parameter_name (@parameter_elements_to_check) {

        ## Test if element from query array exists truth hash
        check_element_exists_in_hash(
            {
                truth_href     => \%parameter,
                queryies       => \@{ $active_parameter{$parameter_name} },
                parameter_name => $parameter_name,
            }
        );
    }

## Parameters with key(s) that have elements as MIP program names
    my @parameter_key_to_check = qw(associated_program);
  PARAMETER:
    foreach my $parameter ( keys %parameter ) {

      KEY:
        foreach my $parameter_name (@parameter_key_to_check) {

            if ( exists $parameter{$parameter}{$parameter_name} ) {

                ## Test if element from query array exists truth hash
                check_element_exists_in_hash(
                    {
                        truth_href => \%parameter,
                        queryies =>
                          \@{ $parameter{$parameter}{$parameter_name} },
                        parameter_name => $parameter_name,
                    }
                );
            }
        }
    }

## Check programs in path, and executable
    check_command_in_path(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
        }
    );

## Test that the family_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
    check_unique_ids(
        {
            active_parameter_href => \%active_parameter,
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        }
    );

## Check sample_id provided in hash parameter is included in the analysis and only represented once
    check_sample_id_in_parameter(
        {
            active_parameter_href => \%active_parameter,
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
            parameter_names_ref =>
              [qw{ analysis_type expected_coverage sample_origin }],
            parameter_href => \%parameter,
        }
    );

## Check sample_id provided in hash path parameter is included in the analysis and only represented once
    check_sample_id_in_parameter_path(
        {
            active_parameter_href => \%active_parameter,
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
            parameter_names_ref   => [qw{ infile_dirs exome_target_bed }],
        }
    );

## Check that VEP directory and VEP cache match
    check_vep_directories(
        {
            vep_directory_path_ref  => \$active_parameter{vep_directory_path},
            vep_directory_cache_ref => \$active_parameter{vep_directory_cache},
        }
    );

## Check that the supplied vcfanno toml frequency file match record 'file=' within toml config file
    if (   ( $active_parameter{psv_combinevariantcallsets} > 0 )
        && ( $active_parameter{sv_vcfanno} > 0 ) )
    {

        check_vcfanno_toml(
            {
                vcfanno_file_toml => $active_parameter{sv_vcfanno_config},
                vcfanno_file_freq => $active_parameter{sv_vcfanno_config_file},
            }
        );
    }

    check_snpsift_keys(
        {
            snpsift_annotation_files_href =>
              \%{ $active_parameter{snpsift_annotation_files} },
            snpsift_annotation_outinfo_key_href =>
              \%{ $active_parameter{snpsift_annotation_outinfo_key} },
        }
    );

## Adds dynamic aggregate information from definitions to parameter hash
    set_dynamic_parameter(
        {
            parameter_href => \%parameter,
            aggregates_ref => [
                ## Collects all programs that MIP can handle
                q{type:program},
                ## Collects all variant_callers
                q{program_type:variant_callers},
                ## Collects all structural variant_callers
                q{program_type:structural_variant_callers},
                ## Collect all aligners
                q{program_type:aligners},
                ## Collects all references in that are supposed to be in reference directory
                q{reference:reference_dir},
            ],
        }
    );

## Check correct value for program mode in MIP
    check_program_mode(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter
        }
    );

## Get initiation program, downstream dependencies and update program modes
    if ( $active_parameter{start_with_program} ) {

        my %dependency_tree = load_yaml(
            {
                yaml_file => catfile(
                    $Bin, qw{ definitions rare_disease_initiation.yaml }
                ),
            }
        );

        my @start_with_programs;
        my $is_program_found = 0;
        my $is_chain_found   = 0;

        ## Collects all downstream programs from initation point
        get_dependency_tree(
            {
                dependency_tree_href => \%dependency_tree,
                is_program_found_ref => \$is_program_found,
                is_chain_found_ref   => \$is_chain_found,
                program              => $active_parameter{start_with_program},
                start_with_programs_ref => \@start_with_programs,
            }
        );

        ## Update program mode depending on start with flag
        update_program_mode_with_start_with(
            {
                active_parameter_href => \%active_parameter,
                programs_ref => \@{ $parameter{dynamic_parameter}{program} },
                start_with_programs_ref => \@start_with_programs,
            }
        );
    }

## Update program mode depending on dry_run_all flag
    update_program_mode_with_dry_run_all(
        {
            active_parameter_href => \%active_parameter,
            programs_ref => \@{ $parameter{dynamic_parameter}{program} },
            dry_run_all  => $active_parameter{dry_run_all},
        }
    );

## Check that the correct number of aligners is used in MIP and sets the aligner flag accordingly
    check_aligner(
        {
            active_parameter_href => \%active_parameter,
            broadcasts_ref        => \@broadcasts,
            parameter_href        => \%parameter,
            verbose               => $active_parameter{verbose},
        }
    );

## Check that all active variant callers have a prioritization order and that the prioritization elements match a supported variant caller.
    my %priority_call_parameter = (
        variant_callers            => 'gatk_combinevariants_prioritize_caller',
        structural_variant_callers => 'sv_svdb_merge_prioritize',
    );
    while ( my ( $variant_caller_type, $prioritize_parameter_name ) =
        each %priority_call_parameter )
    {

        ## Check if we have any active callers
        my $activate_caller_tracker = 0;
        foreach my $variant_caller (
            @{ $parameter{dynamic_parameter}{$variant_caller_type} } )
        {

            if ( $active_parameter{$variant_caller} > 0 ) {

                $activate_caller_tracker++;
            }
        }
        if ( $activate_caller_tracker > 0 ) {

            check_prioritize_variant_callers(
                {
                    parameter_href        => \%parameter,
                    active_parameter_href => \%active_parameter,
                    variant_callers_ref =>
                      \@{ $parameter{dynamic_parameter}{$variant_caller_type} },
                    parameter_names_ref => \$prioritize_parameter_name,
                }
            );
        }
    }

## Broadcast set parameters info
    foreach my $parameter_info (@broadcasts) {

        $log->info($parameter_info);
    }

## Update program mode depending on analysis run value as some programs are not applicable for e.g. wes
    update_program_mode_for_analysis_type(
        {
            active_parameter_href => \%active_parameter,
            consensus_analysis_type =>
              $parameter{dynamic_parameter}{consensus_analysis_type},
            log          => $log,
            programs_ref => [
                qw{ cnvnator delly_call delly_reformat tiddit samtools_subsample_mt }
            ],
        }
    );

## Update prioritize flag depending on analysis run value as some programs are not applicable for e.g. wes
    $active_parameter{sv_svdb_merge_prioritize} = update_prioritize_flag(
        {
            prioritize_key => $active_parameter{sv_svdb_merge_prioritize},
            programs_ref   => [qw{ cnvnator delly_call delly_reformat tiddit }],
            consensus_analysis_type =>
              $parameter{dynamic_parameter}{consensus_analysis_type},
        }
    );

## Write config file for family
    if ( $active_parameter{config_file_analysis} ne 0 ) {

        ## Create directory unless it already exists
        make_path( dirname( $active_parameter{config_file_analysis} ) );

        ## Remove previous analysis specific info not relevant for current run e.g. log file, sample_ids which are read from pedigree or cmd
        my @remove_keys = (qw{ associated_program });

      KEY:
        foreach my $key (@remove_keys) {

            delete $active_parameter{$key};
        }

        ## Writes a YAML hash to file
        write_yaml(
            {
                yaml_href      => \%active_parameter,
                yaml_file_path => $active_parameter{config_file_analysis},
            }
        );
        $log->info( 'Wrote: ' . $active_parameter{config_file_analysis}, "\n" );

        ## Add to qc_sample_info
        $sample_info{config_file_analysis} =
          $active_parameter{config_file_analysis};
    }

## Detect the gender included in current analysis
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
            file_info_href         => \%file_info,
            human_genome_reference => $active_parameter{human_genome_reference},
        }
    );

## Update contigs depending on settings in run (wes or if only male samples)
    update_contigs_for_run(
        {
            file_info_href     => \%file_info,
            analysis_type_href => \%{ $active_parameter{analysis_type} },
            found_male         => $active_parameter{found_male},
        }
    );

## Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
    @{ $file_info{sorted_select_file_contigs} } = size_sort_select_file_contigs(
        {
            file_info_href => \%file_info,
            consensus_analysis_type_ref =>
              \$parameter{dynamic_parameter}{consensus_analysis_type},
            hash_key_to_sort        => 'select_file_contigs',
            hash_key_sort_reference => 'contigs_size_ordered',
        }
    );

    if ( $active_parameter{verbose} ) {
## Write CMD to MIP log file
        write_cmd_mip_log(
            {
                parameter_href        => \%parameter,
                active_parameter_href => \%active_parameter,
                order_parameters_ref  => \@order_parameters,
                script_ref            => \$script,
                log_file_ref          => \$active_parameter{log_file},
                mip_version_ref       => \$VERSION,
            }
        );
    }

## Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any of the files exist
    collect_infiles(
        {
            active_parameter_href => \%active_parameter,
            indir_path_href       => \%indir_path,
            infile_href           => \%infile,
        }
    );

## Reformat files for MIP output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames
    my $uncompressed_file_switch = infiles_reformat(
        {
            active_parameter_href           => \%active_parameter,
            sample_info_href                => \%sample_info,
            file_info_href                  => \%file_info,
            infile_href                     => \%infile,
            indir_path_href                 => \%indir_path,
            infile_lane_prefix_href         => \%infile_lane_prefix,
            infile_both_strands_prefix_href => \%infile_both_strands_prefix,
            lane_href                       => \%lane,
            job_id_href                     => \%job_id,
            outaligner_dir_ref => \$active_parameter{outaligner_dir},
            program_name       => 'infiles_reformat',
        }
    );

## Creates all fileendings as the samples is processed depending on the chain of modules activated
    create_file_endings(
        {
            parameter_href          => \%parameter,
            active_parameter_href   => \%active_parameter,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            order_parameters_ref    => \@order_parameters,
        }
    );

## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            parameter_href        => \%parameter,
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
            execution_mode        => 'system',
            fam_file_path         => catfile(
                $active_parameter{outdata_dir},
                $active_parameter{family_id},
                $active_parameter{family_id} . '.fam'
            ),
        }
    );

## Add to SampleInfo
    add_to_sample_info(
        {
            active_parameter_href => \%active_parameter,
            sample_info_href      => \%sample_info,
            file_info_href        => \%file_info,
        }
    );

############
####MAIN####
############

    if ( not $active_parameter{dry_run_all} ) {

        my %no_dry_run_info = (
            analysisrunstatus => q{not_finished},
            analysis_date     => $date_time_stamp,
            mip_version       => $VERSION,
        );

      KEY_VALUE_PAIR:
        while ( my ( $key, $value ) = each %no_dry_run_info ) {

            $sample_info{$key} = $value;
        }
    }

    my $consensus_analysis_type =
      $parameter{dynamic_parameter}{consensus_analysis_type};

## Split of fastq files in batches
    if ( $active_parameter{psplit_fastq_file} ) {

        $log->info(q{[Split fastq files in batches]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

            ## Split input fastq files into batches of reads, versions and compress. Moves original file to subdirectory
            analysis_split_fastq_file(
                {
                    parameter_href        => \%parameter,
                    active_parameter_href => \%active_parameter,
                    infile_href           => \%infile,
                    job_id_href           => \%job_id,
                    insample_directory    => $indir_path{$sample_id},
                    outsample_directory   => $indir_path{$sample_id},
                    sample_id             => $sample_id,
                    program_name          => q{split_fastq_file},
                    sequence_read_batch =>
                      $active_parameter{split_fastq_file_read_batch},
                }
            );
        }

        ## End here if this module is turned on
        exit;
    }

## GZip of fastq files
    if (   $active_parameter{pgzip_fastq}
        && $uncompressed_file_switch eq q{uncompressed} )
    {

        $log->info(q{[Gzip for fastq files]});

      SAMPLES:
        foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

            ## Determine which sample id had the uncompressed files
          INFILES:
            foreach my $infile ( @{ $infile{$sample_id} } ) {

                my $infile_suffix = $parameter{pgzip_fastq}{infile_suffix};

                if ( $infile =~ /$infile_suffix$/ ) {

                    ## Automatically gzips fastq files
                    analysis_gzip_fastq(
                        {
                            parameter_href          => \%parameter,
                            active_parameter_href   => \%active_parameter,
                            sample_info_href        => \%sample_info,
                            infile_href             => \%infile,
                            infile_lane_prefix_href => \%infile_lane_prefix,
                            job_id_href             => \%job_id,
                            insample_directory      => $indir_path{$sample_id},
                            sample_id               => $sample_id,
                            program_name            => q{gzip_fastq},
                        }
                    );

                    # Call once per sample_id
                    last INFILES;
                }
            }
        }
    }

### Cancer
    if ( $consensus_analysis_type eq q{cancer} )

    {

        $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );

        ## Pipeline recipe for cancer data
        pipeline_cancer(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                indir_path_href         => \%indir_path,
                infile_href             => \%infile,
                infile_lane_prefix_href => \%infile_lane_prefix,
                lane_href               => \%lane,
                job_id_href             => \%job_id,
                outaligner_dir          => $active_parameter{outaligner_dir},
                log                     => $log,
            }
        );
    }

### RNA
    if ( $consensus_analysis_type eq q{wts} ) {

        $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );

        ## Pipeline recipe for rna data
        pipeline_rna(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                indir_path_href         => \%indir_path,
                infile_href             => \%infile,
                infile_lane_prefix_href => \%infile_lane_prefix,
                lane_href               => \%lane,
                job_id_href             => \%job_id,
                outaligner_dir          => $active_parameter{outaligner_dir},
                log                     => $log,
            }
        );
    }

### WES|WGS
    if (   $consensus_analysis_type eq q{wgs}
        || $consensus_analysis_type eq q{wes}
        || $consensus_analysis_type eq q{mixed} )
    {

        $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );

        ## Pipeline recipe for rna data
        pipeline_rare_disease(
            {
                parameter_href          => \%parameter,
                active_parameter_href   => \%active_parameter,
                sample_info_href        => \%sample_info,
                file_info_href          => \%file_info,
                indir_path_href         => \%indir_path,
                infile_href             => \%infile,
                infile_lane_prefix_href => \%infile_lane_prefix,
                lane_href               => \%lane,
                job_id_href             => \%job_id,
                outaligner_dir          => $active_parameter{outaligner_dir},
                log                     => $log,
            }
        );
    }

## Write QC for programs used in analysis
    # Write SampleInfo to yaml file
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

}
######################
####Sub routines######
######################

sub collect_infiles {

##collect_infiles

##Function : Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any files exist.
##Returns  : ""
##Arguments: $active_parameter_href, $indir_path_href, $infile_href
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $indir_path_href       => Indirectories path(s) hash {REF}
##         : $infile_href           => Infiles hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $indir_path_href;
    my $infile_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        indir_path_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$indir_path_href
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    $log->info("Reads from platform:\n");

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
    {    #Collects inputfiles govern by sample_ids

        ## Return the key if the hash value and query match
        my $infile_directory_ref = \get_matching_values_key(
            {
                active_parameter_href => $active_parameter_href,
                query_value_ref       => \$sample_id,
                parameter_name        => "infile_dirs",
            }
        );

        my @infiles;

        ## Collect all fastq files in supplied indirectories
        my $rule = Path::Iterator::Rule->new;
        $rule->skip_subdirs("original_fastq_files")
          ;    #Ignore if original fastq files sub directory
        $rule->name("*.fastq*");    #Only look for fastq or fastq.gz files
        my $it = $rule->iter($$infile_directory_ref);

        while ( my $file = $it->() ) {    #Iterate over directory

            my ( $volume, $directory, $fastq_file ) = splitpath($file);
            push( @infiles, $fastq_file );
        }
        chomp(@infiles);    #Remove newline from every entry in array

        if ( !@infiles ) {  #No "*.fastq*" infiles

            $log->fatal(
"Could not find any '.fastq' files in supplied infiles directory "
                  . $$infile_directory_ref,
                "\n"
            );
            exit 1;
        }
        foreach my $infile (@infiles)
        {    #Check that inFileDirs/infile contains sample_id in filename

            unless ( $infile =~ /$sample_id/ ) {

                $log->fatal(
                    "Could not detect sample_id: "
                      . $sample_id
                      . " in supplied infile: "
                      . $$infile_directory_ref . "/"
                      . $infile,
                    "\n"
                );
                $log->fatal(
"Check that: '--sample_ids' and '--inFileDirs' contain the same sample_id and that the filename of the infile contains the sample_id.",
                    "\n"
                );
                exit 1;
            }
        }
        $log->info( "Sample id: " . $sample_id . "\n" );
        $log->info("\tInputfiles:\n");

        ## Log each file from platform
        foreach my $file (@infiles) {

            $log->info( "\t\t", $file, "\n" );    #Indent for visability
        }
        $indir_path_href->{$sample_id} =
          $$infile_directory_ref;                 #Catch inputdir path
        $infile_href->{$sample_id} = [@infiles];  #Reload files into hash
    }
}

sub infiles_reformat {

## Function : Reformat files for MIP output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames.
## Returns  : "$uncompressed_file_counter"
## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $file_info_href                  => File info hash {REF}
##          : $indir_path_href                 => Indirectories path(s) hash {REF}
##          : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##          : $infile_href                     => Infiles hash {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $job_id_href                     => Job id hash {REF}
##          : $lane_href                       => The lane info hash {REF}
##          : $outaligner_dir_ref              => Outaligner_dir used in the analysis {REF}
##          : $program_name                    => Program name {REF}
##          : $sample_info_href                => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $indir_path_href;
    my $infile_both_strands_prefix_href;
    my $infile_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $lane_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $outaligner_dir_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        indir_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$indir_path_href,
            strict_type => 1,
        },
        infile_both_strands_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_both_strands_prefix_href,
            strict_type => 1,
        },
        infile_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        lane_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$lane_href,
            strict_type => 1,
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir_ref,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
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

    use MIP::Check::Parameter qw{ check_gzipped };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

# Used to decide later if any inputfiles needs to be compressed before starting analysis
    my $uncompressed_file_counter = 0;

  SAMPLE_ID:
    for my $sample_id ( keys %{$infile_href} ) {

        # Needed to be able to track when lanes are finished
        my $lane_tracker = 0;

      INFILE:
        while ( my ( $file_index, $file_name ) =
            each( @{ $infile_href->{$sample_id} } ) )
        {

            ## Check if a file is gzipped.
            my $compressed_switch =
              check_gzipped( { file_name => $file_name, } );
            my $read_file_command = q{zcat};

            ## Not compressed
            if ( not $compressed_switch ) {

                ## File needs compression before starting analysis. Note: All files are rechecked downstream and uncompressed ones are gzipped automatically
                $uncompressed_file_counter = q{uncompressed};
                $read_file_command         = q{cat};
            }

            ## Parse 'new' no "index" format $1=lane, $2=date,
            ## $3=Flow-cell, $4=Sample_id, $5=index,$6=direction
            if (
                $file_name =~ /(\d+)_(\d+)_([^_]+)_([^_]+)_([^_]+)_(\d).fastq/ )
            {

                ## Check that the sample_id provided and sample_id in infile name match.
                check_sample_id_match(
                    {
                        active_parameter_href => $active_parameter_href,
                        file_index            => $file_index,
                        infile_href           => $infile_href,
                        infile_sample_id => $4,    #$4 = Sample_id from filename
                        sample_id => $sample_id,
                    }
                );

                ## Adds information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elementents.
                add_infile_info(
                    {
                        active_parameter_href => $active_parameter_href,
                        compressed_switch     => $compressed_switch,
                        date                  => $2,
                        direction             => $6,
                        file_info_href        => $file_info_href,
                        file_index            => $file_index,
                        flowcell              => $3,
                        index                 => $5,
                        indir_path_href       => $indir_path_href,
                        infile_both_strands_prefix_href =>
                          $infile_both_strands_prefix_href,
                        infile_href             => $infile_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        lane                    => $1,
                        lane_href               => $lane_href,
                        lane_tracker_ref        => \$lane_tracker,
                        sample_id               => $4,
                        sample_info_href        => $sample_info_href,
                    }
                );
            }
            else
            {    #No regexp match i.e. file does not follow filename convention

                $log->warn(
                        q{Could not detect MIP file name convention for file: }
                      . $file_name
                      . q{.} );
                $log->warn(
                    q{Will try to find mandatory information in fastq header.});

                ## Check that file name at least contains sample_id
                if ( $file_name !~ /$sample_id/ ) {

                    $log->fatal(
q{Please check that the file name contains the sample_id.}
                    );
                }

                ## Get run info from fastq file header
                my @fastq_info_headers = get_run_info(
                    {
                        directory         => $indir_path_href->{$sample_id},
                        file              => $file_name,
                        read_file_command => $read_file_command,
                    }
                );

                ## Adds information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elementents.
                add_infile_info(
                    {
                        active_parameter_href => $active_parameter_href,
                        compressed_switch     => $compressed_switch,
                        ## fastq format does not contain a date of the run,
                        ## so fake it with constant impossible date
                        date            => q{000101},
                        direction       => $fastq_info_headers[4],
                        file_index      => $file_index,
                        file_info_href  => $file_info_href,
                        flowcell        => $fastq_info_headers[2],
                        index           => $fastq_info_headers[5],
                        indir_path_href => $indir_path_href,
                        infile_both_strands_prefix_href =>
                          $infile_both_strands_prefix_href,
                        infile_href             => $infile_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        lane                    => $fastq_info_headers[3],
                        lane_href               => $lane_href,
                        lane_tracker_ref        => \$lane_tracker,
                        sample_id               => $sample_id,
                        sample_info_href        => $sample_info_href,
                    }
                );

                $log->info(
                        q{Found following information from fastq header: lane=}
                      . $fastq_info_headers[3]
                      . q{ flow-cell=}
                      . $fastq_info_headers[2]
                      . q{ index=}
                      . $fastq_info_headers[5]
                      . q{ direction=}
                      . $fastq_info_headers[4],
                );
                $log->warn(
q{Will add fake date '20010101' to follow file convention since this is not recorded in fastq header}
                );
            }
        }
    }
    return $uncompressed_file_counter;
}

sub check_sample_id_match {

##check_sample_id_match

##Function : Check that the sample_id provided and sample_id in infile name match.
##Returns  : ""
##Arguments: $active_parameter_href, $infile_href, $sample_id, $infile_sample_id, $file_index
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $infile_href           => Infiles hash {REF}
##         : $sample_id             => Sample id from user
##         : $infile_sample_id      => Sample_id collect with regexp from infile
##         : $file_index            => Counts the number of infiles

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $infile_href;
    my $sample_id;
    my $infile_sample_id;
    my $file_index;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        infile_sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_sample_id
        },
        file_index => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_index
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %seen = ( $infile_sample_id => 1 );    #Add input as first increment

    foreach my $sample_id_supplied ( @{ $active_parameter_href->{sample_ids} } )
    {

        $seen{$sample_id_supplied}++;
    }
    unless ( $seen{$infile_sample_id} > 1 ) {

        $log->fatal( $sample_id
              . " supplied and sample_id "
              . $infile_sample_id
              . " found in file : "
              . $infile_href->{$sample_id}[$file_index]
              . " does not match. Please rename file to match sample_id: "
              . $sample_id
              . "\n" );
        exit 1;
    }
}

sub get_run_info {

##get_run_info

##Function : Get run info from fastq file header
##Returns  : ""
##Arguments: $directory, $read_file, $file
##         : $directory       => Directory of file
##         : $read_file_command => Command used to read file
##         : $file            => File to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $read_file_command;
    my $file;

    my $tmpl = {
        directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$directory
        },
        read_file_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$read_file_command
        },
        file =>
          { required => 1, defined => 1, strict_type => 1, store => \$file },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $fastq_header_regexp =
q?perl -nae 'chomp($_); if($_=~/^(@\w+):(\w+):(\w+):(\w+)\S+\s(\w+):\w+:\w+:(\w+)/) {print $1." ".$2." ".$3." ".$4." ".$5." ".$6."\n";} if($.=1) {last;}' ?;

    my $pwd = cwd();      #Save current direcory
    chdir($directory);    #Move to sample_id infile directory

    my $fastq_info_headers = `$read_file_command $file | $fastq_header_regexp;`
      ;                   #Collect fastq header info
    my @fastq_info_headers = split( " ", $fastq_info_headers );

    chdir($pwd);          #Move back to original directory

    unless ( scalar(@fastq_info_headers) eq 6 ) {

        $log->fatal(
"Could not detect reuired sample sequencing run info from fastq file header - PLease proved MIP file in MIP file convention format to proceed\n"
        );
        exit 1;
    }

    return @fastq_info_headers;
}

sub add_infile_info {

##add_infile_info

##Function : Adds information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elementents.
##Returns  : ""
##Arguments: $active_parameter_href, $sample_info_href, $file_info_href, $infile_href, $infile_lane_prefix_href, $infile_both_strands_prefix_href, $indir_path_href, $lane_href, $lane, $date, $flowcell, $sample_id, $index, $direction, $lane_tracker_ref, $file_index, $compressed_switch
##         : $active_parameter_href              => Active parameters for this analysis hash {REF}
##         : $sample_info_href                   => Info on samples and family hash {REF}
##         : $file_info_href                     => File info hash {REF}
##         : $infile_href                        => Infiles hash {REF}
##         : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##         : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##         : $indir_path_href                    => Indirectories path(s) hash {REF}
##         : $lane_href                          => The lane info hash {REF}
##         : $lane                               => Flow-cell lane
##         : $date                               => Flow-cell sequencing date
##         : $flowcell                           => Flow-cell id
##         : $sample_id                          => Sample id
##         : $index                              => The DNA library preparation molecular barcode
##         : $direction                          => Sequencing read direction
##         : $lane_tracker_ref                   => Counts the number of lanes sequenced {REF}
##         : $file_index                         => Index of file
##         : $compressed_switch                  => ".fastq.gz" or ".fastq" info governs zcat or cat downstream

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_href;
    my $indir_path_href;
    my $infile_lane_prefix_href;
    my $infile_both_strands_prefix_href;
    my $lane_href;
    my $lane_tracker_ref;
    my $sample_id;
    my $lane;
    my $date;
    my $flowcell;
    my $index;
    my $direction;
    my $file_index;
    my $compressed_switch;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_href
        },
        indir_path_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$indir_path_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        infile_both_strands_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_both_strands_prefix_href
        },
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        lane => {
            required    => 1,
            defined     => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$lane
        },
        lane_tracker_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$lane_tracker_ref
        },
        date =>
          { required => 1, defined => 1, strict_type => 1, store => \$date },
        flowcell => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$flowcell
        },
        index =>
          { required => 1, defined => 1, strict_type => 1, store => \$index },
        direction => {
            required    => 1,
            defined     => 1,
            allow       => [ 1, 2 ],
            strict_type => 1,
            store       => \$direction
        },
        file_index => {
            required    => 1,
            defined     => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$file_index
        },
        compressed_switch => {
            required    => 1,
            defined     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$compressed_switch
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $read_file;
    my $file_at_lane_level_ref;
    my $file_at_direction_level_ref;

    my $parsed_date = Time::Piece->strptime( $date, "%y%m%d" );
    $parsed_date = $parsed_date->ymd;

    if ($compressed_switch) {

        $read_file = "zcat";    #Read file in compressed format
    }
    else {

        $read_file = "cat";     #Read file in uncompressed format
    }

    if ( $direction == 1 ) {    #Read 1

        push( @{ $lane_href->{$sample_id} }, $lane );    #Lane
        $infile_lane_prefix_href->{$sample_id}[$$lane_tracker_ref] =
            $sample_id . "."
          . $date . "_"
          . $flowcell . "_"
          . $index . ".lane"
          . $lane
          ; #Save new format (sample_id_date_flow-cell_index_lane) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).

        $file_at_lane_level_ref =
          \$infile_lane_prefix_href->{$sample_id}[$$lane_tracker_ref];    #Alias
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {sequence_run_type} = "single_end"; #Single_end until proven otherwise

        ## Collect read length from an infile
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {sequence_length} = collect_read_length(
            {
                directory         => $indir_path_href->{$sample_id},
                read_file_command => $read_file,
                file              => $infile_href->{$sample_id}[$file_index],
            }
          );

        ## Check if fastq file is interleaved
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {interleaved} = detect_interleaved(
            {
                directory         => $indir_path_href->{$sample_id},
                read_file_command => $read_file,
                file              => $infile_href->{$sample_id}[$file_index],
            }
          );

        ## Detect "regexp" in string
        $file_info_href->{undetermined_in_file_name}
          { $infile_lane_prefix_href->{$sample_id}[$$lane_tracker_ref] } =
          check_string(
            {
                string => $flowcell,
                regexp => "Undetermined",
            }
          );
        $$lane_tracker_ref++;
    }
    if ( $direction == 2 ) {    #2nd read direction

        $file_at_lane_level_ref =
          \$infile_lane_prefix_href->{$sample_id}[ $$lane_tracker_ref - 1 ]
          ;                     #Alias
        $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
          {sequence_run_type} = 'paired-end'
          ;    #$lane_tracker -1 since it gets incremented after direction eq 1.
    }

    $infile_both_strands_prefix_href->{$sample_id}[$file_index] =
        $sample_id . "."
      . $date . "_"
      . $flowcell . "_"
      . $index . ".lane"
      . $lane . "_"
      . $direction
      ; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).

    $file_at_direction_level_ref =
      \$infile_both_strands_prefix_href->{$sample_id}[$file_index];    #Alias
    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{original_file_name}
      = $infile_href->{$sample_id}[$file_index];    #Original file_name

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}
      {original_file_name_prefix} =
        $lane . "_"
      . $date . "_"
      . $flowcell . "_"
      . $sample_id . "_"
      . $index . "_"
      . $direction;    #Original file_name, but no ending

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{lane} =
      $lane;           #Save sample lane

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{date} =
      $parsed_date;    #Save Sequence run date

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{flowcell} =
      $flowcell;       #Save Sequence flow-cell

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{sample_barcode} =
      $index;          #Save sample barcode

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{run_barcode} =
      $date . "_" . $flowcell . "_" . $lane . "_" . $index;    #Save run barcode

    $sample_info_href->{sample}{$sample_id}{file}{$$file_at_lane_level_ref}
      {read_direction_file}{$$file_at_direction_level_ref}{read_direction} =
      $direction;
}

sub detect_interleaved {

##detect_interleaved

##Function : Detect if fastq file is interleaved
##Returns  : "1(=interleaved)"
##Arguments: $directory, $read_file, $file
##         : $directory         => Directory of file
##         : $read_file_command => Command used to read file
##         : $file              => File to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $read_file_command;
    my $file;

    my $tmpl = {
        directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$directory
        },
        read_file_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$read_file_command
        },
        file =>
          { required => 1, defined => 1, strict_type => 1, store => \$file },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $interleaved_regexp =
q?perl -nae 'chomp($_); if( ($_=~/^@\S+:\w+:\w+:\w+\S+\s(\w+):\w+:\w+:\w+/) && ($.==5) ) {print $1."\n";last;} elsif ($.==6) {last;}' ?;

    my $pwd = cwd();      #Save current direcory
    chdir($directory);    #Move to sample_id infile directory

    my $fastq_info_headers = `$read_file_command $file | $interleaved_regexp;`
      ;                   #Collect interleaved info

    if ( !$fastq_info_headers ) {

        my $interleaved_regexp =
q?perl -nae 'chomp($_); if( ($_=~/^@\w+-\w+:\w+:\w+:\w+:\w+:\w+:\w+\/(\w+)/) && ($.==5) ) {print $1."\n";last;} elsif ($.==6) {last;}' ?;
        $fastq_info_headers = `$read_file_command $file | $interleaved_regexp;`
          ;               #Collect interleaved info
    }

    chdir($pwd);          #Move back to original directory

    unless ( $fastq_info_headers =~ /[1, 2, 3]/ ) {

        $log->fatal("Malformed fastq file!\n");
        $log->fatal( "Read direction is: "
              . $fastq_info_headers
              . " allowed entries are '1', '2', '3'. Please check fastq file\n"
        );
        exit 1;
    }
    if ( $fastq_info_headers > 1 ) {

        $log->info( "Found interleaved fastq file: " . $file, "\n" );
        return 1;
    }
    return;
}

sub create_file_endings {

## Function : Creates the file_tags depending on which modules are used by the user to relevant chain.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id_ref           => Family id {REF}
##          : $file_info_href          => Info on files hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $order_parameters_ref    => Order of addition to parameter array {REF}
##          : $parameter_href          => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $order_parameters_ref;
    my $parameter_href;

    ## Default(s)
    my $family_id_ref;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            store       => \$family_id_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};

    ## Used to enable seqential build-up of file_tags between modules
    my %temp_file_ending;

  PARAMETER:
    foreach my $order_parameter_element (@$order_parameters_ref) {

        ## Only active parameters
        if ( defined $active_parameter_href->{$order_parameter_element} ) {

            ## Only process programs
            if (
                any { $_ eq $order_parameter_element }
                @{ $parameter_href->{dynamic_parameter}{program} }
              )
            {

                ## MAIN chain
                if ( $parameter_href->{$order_parameter_element}{chain} eq
                    q{MAIN} )
                {

                    ##  File_tag exist
                    if ( $parameter_href->{$order_parameter_element}{file_tag}
                        ne q{nofile_tag} )
                    {

                        ## Alias
                        my $file_ending_ref =
                          \$parameter_href->{$order_parameter_element}
                          {file_tag};

                        ###MAIN/Per sample_id
                      SAMPLE_ID:
                        foreach my $sample_id (
                            @{ $active_parameter_href->{sample_ids} } )
                        {

                            ## File_ending should be added
                            if ( $active_parameter_href
                                ->{$order_parameter_element} > 0 )
                            {

                                ## Special case
                                if ( $order_parameter_element eq
                                    q{ppicardtools_mergesamfiles} )
                                {

                                    $file_info_href->{$sample_id}
                                      {ppicardtools_mergesamfiles}{file_tag} =
                                      $temp_file_ending{$sample_id} . "";
                                }
                                else {

                                    if ( defined $temp_file_ending{$sample_id} )
                                    {

                                        $file_info_href->{$sample_id}
                                          {$order_parameter_element}{file_tag}
                                          = $temp_file_ending{$sample_id}
                                          . $$file_ending_ref;
                                    }
                                    else {
                                        ## First module that should add filending

                                        $file_info_href->{$sample_id}
                                          {$order_parameter_element}{file_tag}
                                          = $$file_ending_ref;
                                    }
                                }
                            }
                            else {
                                ## Do not add new module file_tag

                                $file_info_href->{$sample_id}
                                  {$order_parameter_element}{file_tag} =
                                  $temp_file_ending{$sample_id};
                            }

                            ## To enable sequential build-up of fileending
                            $temp_file_ending{$sample_id} =
                              $file_info_href->{$sample_id}
                              {$order_parameter_element}{file_tag};
                        }

                        ###MAIN/Per family_id
                        ## File_ending should be added
                        if ( $active_parameter_href->{$order_parameter_element}
                            > 0 )
                        {

                            ## Special case - do nothing
                            if ( $order_parameter_element eq
                                q{ppicardtools_mergesamfiles} )
                            {
                            }
                            else {

                                if (
                                    defined $temp_file_ending{$$family_id_ref} )
                                {

                                    $file_info_href->{$$family_id_ref}
                                      {$order_parameter_element}{file_tag} =
                                        $temp_file_ending{$$family_id_ref}
                                      . $$file_ending_ref;
                                }
                                else {
                                    ## First module that should add filending

                                    $file_info_href->{$$family_id_ref}
                                      {$order_parameter_element}{file_tag} =
                                      $$file_ending_ref;
                                }

                                ## To enable sequential build-up of fileending
                                $temp_file_ending{$$family_id_ref} =
                                  $file_info_href->{$$family_id_ref}
                                  {$order_parameter_element}{file_tag};
                            }
                        }
                        else {
                            ## Do not add new module file_tag

                            $file_info_href->{$$family_id_ref}
                              {$order_parameter_element}{file_tag} =
                              $temp_file_ending{$$family_id_ref};
                        }
                    }
                }

                ## Other chain(s)
                if ( $parameter_href->{$order_parameter_element}{chain} ne
                    q{MAIN} )
                {

                    ## Alias
                    my $chain_fork =
                      $parameter_href->{$order_parameter_element}{chain};

                    ## File_tag exist
                    if ( $parameter_href->{$order_parameter_element}{file_tag}
                        ne q{nofile_tag} )
                    {

                        ## Alias
                        my $file_ending_ref =
                          \$parameter_href->{$order_parameter_element}
                          {file_tag};

                        ###OTHER/Per sample_id
                      SAMPLE_ID:
                        foreach my $sample_id (
                            @{ $active_parameter_href->{sample_ids} } )
                        {

                            ## File_ending should be added
                            if ( $active_parameter_href
                                ->{$order_parameter_element} > 0 )
                            {

                                if (
                                    not
                                    defined $temp_file_ending{$chain_fork}
                                    {$sample_id} )
                                {

                                    ## Inherit current MAIN chain.
                                    $temp_file_ending{$chain_fork}{$sample_id}
                                      = $temp_file_ending{$sample_id};
                                }
                                if (
                                    defined $temp_file_ending{$chain_fork}
                                    {$sample_id} )
                                {

                                    $file_info_href->{$sample_id}
                                      {$order_parameter_element}{file_tag} =
                                      $temp_file_ending{$chain_fork}{$sample_id}
                                      . $$file_ending_ref;
                                }
                                else {
                                    ## First module that should add filending

                                    $file_info_href->{$sample_id}
                                      {$order_parameter_element}{file_tag} =
                                      $$file_ending_ref;
                                }
                            }
                            else {
                                ## Do not add new module file_tag

                                $file_info_href->{$sample_id}
                                  {$order_parameter_element}{file_tag} =
                                  $temp_file_ending{$chain_fork}{$sample_id};
                            }

                            ## To enable sequential build-up of fileending
                            $temp_file_ending{$chain_fork}{$sample_id} =
                              $file_info_href->{$sample_id}
                              {$order_parameter_element}{file_tag};
                        }
                        ###Other/Per family_id

                        ## File ending should be added
                        if ( $active_parameter_href->{$order_parameter_element}
                            > 0 )
                        {

                            if (
                                not defined $temp_file_ending{$chain_fork}
                                {$$family_id_ref} )
                            {

                                ## Inherit current MAIN chain.
                                $temp_file_ending{$chain_fork}{$$family_id_ref}
                                  = $temp_file_ending{$$family_id_ref};
                            }
                            if (
                                defined $temp_file_ending{$chain_fork}
                                {$$family_id_ref} )
                            {

                                $file_info_href->{$$family_id_ref}
                                  {$order_parameter_element}{file_tag} =
                                  $temp_file_ending{$chain_fork}
                                  {$$family_id_ref} . $$file_ending_ref;
                            }
                            else {
                                ## First module that should add filending

                                $file_info_href->{$$family_id_ref}
                                  {$order_parameter_element}{file_tag} =
                                  $$file_ending_ref;
                            }

                            ## To enable sequential build-up of fileending
                            $temp_file_ending{$chain_fork}{$$family_id_ref} =
                              $file_info_href->{$$family_id_ref}
                              {$order_parameter_element}{file_tag};
                        }
                        else {
                            ## Do not add new module file_tag

                            $file_info_href->{$$family_id_ref}
                              {$order_parameter_element}{file_tag} =
                              $temp_file_ending{$chain_fork}{$$family_id_ref};
                        }
                    }
                }
            }
        }
    }
    return;
}

sub write_cmd_mip_log {

##write_cmd_mip_log

##Function : Write CMD to MIP log file
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $order_parameters_ref, $script_ref, $log_file_ref
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $order_parameters_ref  => Order of addition to parameter array {REF}
##         : $script_ref            => The script that is being executed {REF}
##         : $log_file_ref          => The log file {REF}
##         : $mip_version_ref       => The MIP version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $order_parameters_ref;
    my $script_ref;
    my $log_file_ref;
    my $mip_version_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        order_parameters_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$order_parameters_ref
        },
        script_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$script_ref
        },
        log_file_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$log_file_ref
        },
        mip_version_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$mip_version_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $cmd_line = $$script_ref . " ";

    my @nowrite = (
        "mip",                  "bwa_build_reference",
        "pbamcalibrationblock", "pvariantannotationblock",
        q{associated_program},  q{rtg_build_reference},
    );

  PARAMETER_KEY:
    foreach my $order_parameter_element ( @{$order_parameters_ref} ) {

        if ( defined $active_parameter_href->{$order_parameter_element} ) {

            ## If no config file do not print
            if (   $order_parameter_element eq q{config_file}
                && $active_parameter_href->{config_file} eq 0 )
            {
            }
            else {

                ## If element is part of array - do nothing
                if ( any { $_ eq $order_parameter_element } @nowrite ) {
                }
                elsif (
                    ## Array reference
                    (
                        exists $parameter_href->{$order_parameter_element}
                        {data_type}
                    )
                    && ( $parameter_href->{$order_parameter_element}{data_type}
                        eq q{ARRAY} )
                  )
                {

                    my $separator = $parameter_href->{$order_parameter_element}
                      {element_separator};
                    $cmd_line .= "-"
                      . $order_parameter_element . " "
                      . join(
                        $separator,
                        @{
                            $active_parameter_href->{$order_parameter_element}
                        }
                      ) . " ";
                }
                elsif (
                    ## HASH reference
                    (
                        exists $parameter_href->{$order_parameter_element}
                        {data_type}
                    )
                    && ( $parameter_href->{$order_parameter_element}{data_type}
                        eq q{HASH} )
                  )
                {

                    # First key
                    $cmd_line .= "-" . $order_parameter_element . " ";
                    $cmd_line .= join(
                        "-" . $order_parameter_element . " ",
                        map {
"$_=$active_parameter_href->{$order_parameter_element}{$_} "
                        } (
                            keys %{
                                $active_parameter_href
                                  ->{$order_parameter_element}
                            }
                        )
                    );
                }
                else {

                    $cmd_line .= "-"
                      . $order_parameter_element . " "
                      . $active_parameter_href->{$order_parameter_element}
                      . " ";
                }
            }
        }
    }
    $log->info( $cmd_line,                            "\n" );
    $log->info( q{MIP Version: } . $$mip_version_ref, "\n" );
    $log->info(
        q{Script parameters and info from }
          . $$script_ref
          . q{ are saved in file: }
          . $$log_file_ref,
        "\n"
    );
    return;
}

sub check_unique_ids {

##check_unique_ids

##Function : Test that the family_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
##Returns  : ""
##Arguments: $active_parameter_href, $sample_ids_ref
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_ids_ref        => Array to loop in for parameter {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %seen;    #Hash to test duplicate sample_ids later

    if ( !@$sample_ids_ref ) {

        $log->fatal("Please provide sample_id(s)\n");
        exit 1;
    }

    foreach my $sample_id ( @{$sample_ids_ref} ) {

        $seen{$sample_id}++;    #Increment instance to check duplicates later

        if ( $$family_id_ref eq $sample_id )
        {                       #Family_id cannot be the same as sample_id

            $log->fatal( "Family_id: "
                  . $$family_id_ref
                  . " equals sample_id: "
                  . $sample_id
                  . ". Please make sure that the family_id and sample_id(s) are unique.\n"
            );
            exit 1;
        }
        if ( $seen{$sample_id} > 1 ) {    #Check sample_id are unique

            $log->fatal( "Sample_id: " . $sample_id . " is not uniqe.\n" );
            exit 1;
        }
        if ( $sample_id =~ /_/ )
        { #Sample_id contains "_", which is not allowed according to filename conventions

            $log->fatal( "Sample_id: "
                  . $sample_id
                  . " contains '_'. Please rename sample_id according to MIP's filename convention, removing the '_'.\n"
            );
            exit 1;
        }
    }
}

sub size_sort_select_file_contigs {

## Function : Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
## Returns  : @sorted_contigs
## Arguments: $consensus_analysis_type_ref => Consensus analysis_type {REF}
##          : $file_info_href              => File info hash {REF}
##          : $hash_key_sort_reference     => The hash keys sort reference
##          : $hash_key_to_sort            => The keys to sort

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type_ref;
    my $file_info_href;
    my $hash_key_sort_reference;
    my $hash_key_to_sort;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        consensus_analysis_type_ref => {
            default     => \$$,
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type_ref,
            strict_type => 1,
        },
        hash_key_to_sort => {
            defined     => 1,
            required    => 1,
            store       => \$hash_key_to_sort,
            strict_type => 1,
        },
        hash_key_sort_reference => {
            defined     => 1,
            required    => 1,
            store       => \$hash_key_sort_reference,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Hash qw{ check_element_exist_hash_of_array };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @sorted_contigs;

    ## Sort the contigs depending on reference array
    if ( $file_info_href->{$hash_key_to_sort} ) {

        foreach my $element ( @{ $file_info_href->{$hash_key_sort_reference} } )
        {

            if (
                not check_element_exist_hash_of_array(
                    {
                        element  => $element,
                        hash_ref => $file_info_href,
                        key      => $hash_key_to_sort,
                    }
                )
              )
            {

                push @sorted_contigs, $element;
            }
        }
    }

    ## Test if all contigs collected from select file was sorted by reference contig array
    if ( @sorted_contigs
        && scalar @{ $file_info_href->{$hash_key_to_sort} } !=
        scalar @sorted_contigs )
    {

        foreach my $element ( @{ $file_info_href->{$hash_key_to_sort} } ) {

            ## If element is not part of array
            if ( not any { $_ eq $element } @sorted_contigs ) {

                ## Special case when analysing wes since Mitochondrial contigs have no baits in exome capture kits
                unless ( $$consensus_analysis_type_ref eq q{wes}
                    && $element =~ /MT$|M$/ )
                {

                    $log->fatal( q{Could not detect '##contig'= }
                          . $element
                          . q{ from meta data header in '-vcfparser_select_file' in reference contigs collected from '-human_genome_reference'}
                    );
                    exit 1;
                }
            }
        }
    }
    return @sorted_contigs;
}

sub detect_sample_id_gender {

##detect_sample_id_gender

##Function : Detect gender of the current analysis
##Returns  : "$found_male $found_female $found_other"
##Arguments: $active_parameter_href, $sample_info_href
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $found_male        = 0;
    my $found_female      = 0;
    my $found_other       = 0;
    my $found_other_count = 0;

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        if ( $sample_info_href->{sample}{$sample_id}{sex} =~ /1|^male/ ) { #Male

            $found_male = 1;                                               #Male
        }
        elsif ( $sample_info_href->{sample}{$sample_id}{sex} =~ /2|female/ )
        {    #Female

            $found_female = 1;
        }
        else {    #Other

            $found_male =
              1;    #Include since it might be male to enable analysis of Y.
            $found_other = 1;
            $found_other_count++;
        }
    }
    return $found_male, $found_female, $found_other, $found_other_count;
}

sub check_command_in_path {

##check_command_in_path

##Function : Checking commands in your path and executable
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Unix qw{check_binary_in_path};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    # Track program paths that have already been checked
    my %seen;

    foreach my $parameter_name ( keys %$active_parameter_href ) {

        if (   ( exists( $parameter_href->{$parameter_name}{type} ) )
            && ( $parameter_href->{$parameter_name}{type} eq "program" ) )
        {

            my $program_name_paths_ref =
              \@{ $parameter_href->{$parameter_name}{program_name_path} }
              ;    #Alias

            if (   (@$program_name_paths_ref)
                && ( $active_parameter_href->{$parameter_name} > 0 ) )
            {      #Only check path(s) for active programs

                foreach my $program ( @{$program_name_paths_ref} ) {

                    unless ( $seen{$program} ) {

                        $seen{$program} = check_binary_in_path(
                            {
                                binary => $program,
                                log    => $log,
                            }
                        );
                    }
                }
            }
        }
    }
}

sub add_to_sample_info {

##add_to_sample_info

##Function : Adds parameter info to sample_info
##Returns  : ""
##Arguments: $active_parameter_href, $sample_info_href, $file_info_href, $family_id_ref
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}
##         : $file_info_href        => File info hash {REF}
##         : $family_id_ref         => The family_id_ref {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $human_genome_reference_ref;
    my $outdata_dir;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref,
        },
        human_genome_reference_ref => {
            default =>
              \$arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference_ref
        },
        outdata_dir => {
            default     => $arg_href->{active_parameter_href}{outdata_dir},
            strict_type => 1,
            store       => \$outdata_dir
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::QC::Record qw(add_program_outfile_to_sample_info);

    if ( exists( $active_parameter_href->{analysis_type} ) ) {

        $sample_info_href->{analysis_type} =
          $active_parameter_href->{analysis_type};
    }
    if ( exists( $active_parameter_href->{expected_coverage} ) ) {

        $sample_info_href->{expected_coverage} =
          $active_parameter_href->{expected_coverage};
    }
    if ( exists $active_parameter_href->{sample_origin} ) {

        $sample_info_href->{sample_origin} =
          $active_parameter_href->{sample_origin};
    }
    if ( exists $active_parameter_href->{gatk_path}
        && $active_parameter_href->{gatk_path} )
    {

        my $gatk_version;
        if ( $active_parameter_href->{gatk_path} =~ /GenomeAnalysisTK-([^,]+)/ )
        {

            $gatk_version = $1;
        }
        else {
            # Fall back on actually calling program

            my $jar_path = catfile( $active_parameter_href->{gatk_path},
                "GenomeAnalysisTK.jar" );
            $gatk_version = (`java -jar $jar_path --version 2>&1`);
            chomp $gatk_version;
        }
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'gatk',
                version          => $gatk_version,
            }
        );
    }
    if ( exists $active_parameter_href->{picardtools_path}
        && $active_parameter_href->{picardtools_path} )
    {
        ## To enable addition of version to sample_info

        my $picardtools_version;
        if ( $active_parameter_href->{picardtools_path} =~
            /picard-tools-([^,]+)/ )
        {

            $picardtools_version = $1;
        }
        else {    #Fall back on actually calling program

            my $jar_path = catfile( $active_parameter_href->{picardtools_path},
                q{picard.jar} );
            $picardtools_version =
              (`java -jar $jar_path CreateSequenceDictionary --version 2>&1`);
            chomp $picardtools_version;
        }

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => 'picardtools',
                version          => $picardtools_version,
            }
        );
    }
    my @sambamba_programs =
      ( "pbwa_mem", "psambamba_depth", "markduplicates_sambamba_markdup" );
    foreach my $program (@sambamba_programs) {

        if (   ( defined $active_parameter_href->{$program} )
            && ( $active_parameter_href->{$program} == 1 ) )
        {

            if ( !$active_parameter_href->{dry_run_all} ) {

                my $regexp =
                  q?perl -nae 'if($_=~/sambamba\s(\S+)/) {print $1;last;}'?;
                my $sambamba_version = (`sambamba 2>&1 | $regexp`);
                chomp $sambamba_version;
                add_program_outfile_to_sample_info(
                    {
                        sample_info_href => $sample_info_href,
                        program_name     => 'sambamba',
                        version          => $sambamba_version,
                    }
                );
                last;    #Only need to check once
            }
        }
    }
    if ( exists $active_parameter_href->{pcnvnator} )
    {                    #To enable addition of version to sample_info

        if (   ( $active_parameter_href->{pcnvnator} == 1 )
            && ( !$active_parameter_href->{dry_run_all} ) )
        {

            my $regexp =
              q?perl -nae 'if($_=~/CNVnator\s+(\S+)/) {print $1;last;}'?;
            my $cnvnator_version = (`cnvnator 2>&1 | $regexp`);
            chomp $cnvnator_version;
            add_program_outfile_to_sample_info(
                {
                    sample_info_href => $sample_info_href,
                    program_name     => 'cnvnator',
                    version          => $cnvnator_version,
                }
            );
        }
    }
    if ( defined($$human_genome_reference_ref) )
    {    #To enable addition of version to sample_info

        $sample_info_href->{human_genome_build}{path} =
          $$human_genome_reference_ref;
        $sample_info_href->{human_genome_build}{source} =
          $file_info_href->{human_genome_reference_source};
        $sample_info_href->{human_genome_build}{version} =
          $file_info_href->{human_genome_reference_version};
    }
    if ( exists( $active_parameter_href->{pedigree_file} ) ) {

        ## Add pedigree_file to sample_info
        $sample_info_href->{pedigree_file}{path} =
          $active_parameter_href->{pedigree_file};
    }
    if ( exists( $active_parameter_href->{log_file} ) ) {

        my $path = dirname( dirname( $active_parameter_href->{log_file} ) );
        $sample_info_href->{log_file_dir} =
          $path;    #Add log_file_dir to SampleInfoFile
        $sample_info_href->{last_log_file_path} =
          $active_parameter_href->{log_file};
    }
}

sub check_vep_directories {

##check_vep_directories

##Function : Compare VEP directory and VEP chache versions
##Returns  : ""
##Arguments: $vep_directory_path_ref, $vep_directory_cache_ref
##         : $vep_directory_path_ref  => VEP directory path {REF}
##         : $vep_directory_cache_ref => VEP cache directory path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_directory_path_ref;
    my $vep_directory_cache_ref;

    my $tmpl = {
        vep_directory_path_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$vep_directory_path_ref
        },
        vep_directory_cache_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$vep_directory_cache_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    if ( $$vep_directory_path_ref =~ /ensembl-tools-release-(\d+)/ ) {

        my $vep_directory_path_version = $1;

        unless ( $$vep_directory_cache_ref =~
            /ensembl-tools-release-$vep_directory_path_version/ )
        {

            print $log->fatal(
                "Differing versions between '-vep_directory_path': "
                  . $$vep_directory_path_ref
                  . " and '-vep_directory_cache': "
                  . $$vep_directory_cache_ref,
                "\n"
            );
            exit 1;
        }
    }

}

sub detect_founders {

##detect_founders

##Function : Detect number of founders (i.e. parents ) based on pedigree file
##Returns  : ""|1
##Arguments: $active_parameter_href,
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @founders;

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        my $father_info =
          $sample_info_href->{sample}{$sample_id}{father};    #Alias
        my $mother_info =
          $sample_info_href->{sample}{$sample_id}{mother};    #Alias

        if ( ( defined($father_info) ) && ( $father_info ne 0 ) ) {    #Child

            if (
                any { $_ eq $father_info }
                @{ $active_parameter_href->{sample_ids} }
              )
            {    #If element is part of array

                push( @founders, $father_info );
            }
        }
        if ( ( defined($mother_info) ) && ( $mother_info ne 0 ) ) {    #Child

            if (
                any { $_ eq $mother_info }
                @{ $active_parameter_href->{sample_ids} }
              )
            {    #If element is part of array

                push( @founders, $mother_info );
            }
        }
    }
    return scalar(@founders);
}

sub detect_trio {

##detect_trio

##Function : Detect family constellation based on pedigree file
##Returns  : ""|1
##Arguments: $active_parameter_href, $sample_info_href
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %trio;

    if ( scalar( @{ $active_parameter_href->{sample_ids} } ) eq 1 ) {

        $log->info(
            "Found single sample: " . $active_parameter_href->{sample_ids}[0],
            "\n" );
        return;
    }
    elsif ( scalar( @{ $active_parameter_href->{sample_ids} } ) eq 3 ) {

        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $father_info =
              $sample_info_href->{sample}{$sample_id}{father};    #Alias
            my $mother_info =
              $sample_info_href->{sample}{$sample_id}{mother};    #Alias

            if ( ( $father_info ne 0 ) && ( $mother_info ne 0 ) ) {    #Child

                $trio{child} = $sample_id;

                if (
                    any { $_ eq $father_info }
                    @{ $active_parameter_href->{sample_ids} }
                  )
                {    #If element is part of array

                    $trio{father} = $father_info;
                }
                if (
                    any { $_ eq $mother_info }
                    @{ $active_parameter_href->{sample_ids} }
                  )
                {    #If element is part of array

                    $trio{mother} = $mother_info;
                }
            }
        }
        if ( scalar( keys %trio ) == 3 ) {

            $log->info(
                "Found trio: Child = "
                  . $trio{child}
                  . ", Father = "
                  . $trio{father}
                  . ", Mother = "
                  . $trio{mother},
                "\n"
            );
            return 1;
        }
    }
}

sub check_string {

##check_string

##Function : Detect "regexp" in string
##Returns  : ""|1
##Arguments: $string
##         : $string => String to be searched
##         : $regexp => regexp to use on string

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $string;
    my $regexp;

    my $tmpl = {
        string =>
          { required => 1, defined => 1, strict_type => 1, store => \$string },
        regexp =>
          { required => 1, defined => 1, strict_type => 1, store => \$regexp },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $string =~ /$regexp/ ) {

        return 1;
    }
}

sub check_prioritize_variant_callers {

## Function : Check that all active variant callers have a prioritization order and that the prioritization elements match a supported variant caller.
## Returns  :
## Arguments: $parameter_href        => Parameter hash {REF}
##          : $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $variant_callers_ref   => Variant callers to check {REF}
##          : $parameter_names_ref   => Parameter name list {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $variant_callers_ref;
    my $parameter_names_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
            ,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
            ,
        },
        parameter_names_ref => {
            required => 1,
            defined  => 1,
            default  => [],
            store    => \$parameter_names_ref,
        },
        variant_callers_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$variant_callers_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @priority_calls =
      split( ",", $active_parameter_href->{$$parameter_names_ref} );

    ## No matching variant caller
    my @variant_caller_aliases;

    ## Check that all active variant callers have a priority order
  CALLER:
    foreach my $variant_caller ( @{$variant_callers_ref} ) {

        my $variant_caller_alias =
          $parameter_href->{$variant_caller}{outdir_name};
        push @variant_caller_aliases, $variant_caller_alias;

        ## Only active programs
        if ( $active_parameter_href->{$variant_caller} > 0 ) {

            ## If element is not part of string
            if ( !( any { $_ eq $variant_caller_alias } @priority_calls ) ) {

                $log->fatal( $$parameter_names_ref
                      . q{ does not contain active variant caller: '}
                      . $variant_caller_alias
                      . q{'} );
                exit 1;
            }
        }
        ## Only NOT active programs
        if ( $active_parameter_href->{$variant_caller} == 0 ) {

            ## If element is part of string
            if ( ( any { $_ eq $variant_caller_alias } @priority_calls ) ) {

                $log->fatal( $$parameter_names_ref
                      . q{ contains deactivated variant caller: '}
                      . $variant_caller_alias
                      . q{'} );
                exit 1;
            }
        }
    }

    ## Check that prioritize string contains valid variant call names
    foreach my $prioritize_call (@priority_calls) {

        if ( !( any { $_ eq $prioritize_call } @variant_caller_aliases ) )
        {    #If element is not part of string

            $log->fatal( $$parameter_names_ref . ": '"
                  . $prioritize_call
                  . "' does not match any supported variant caller: '"
                  . join( ",", @variant_caller_aliases )
                  . "'" );
            exit 1;
        }
    }
}

sub check_aligner {

##check_aligner

##Function : Check that the correct number of aligners is used in MIP and sets the outaligner_dir flag accordingly.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $broadcasts_ref, $outaligner_dir_ref
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##         : $outaligner_dir_ref    => Outaligner_dir used in the analysis {REF}
##         : $verbose               => Verbosity level

    my ($arg_href) = @_;

    ## Default(s)
    my $outaligner_dir_ref;
    my $verbose;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $broadcasts_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        broadcasts_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$broadcasts_ref
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref,
        },
        verbose => {
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %aligner;

    foreach my $aligner ( @{ $parameter_href->{dynamic_parameter}{aligners} } )
    {

        if ( $active_parameter_href->{$aligner} > 0 ) {    #Active aligner

            $aligner{total_active_aligner_count}++;
            push( @{ $aligner{active_aligners} }, $aligner );
            $parameter_href->{active_aligner} =
              $aligner;    #Save the active aligner for downstream use

            if ( not defined $$outaligner_dir_ref ) {

                $$outaligner_dir_ref = $parameter_href->{$aligner}{outdir_name}
                  ;    #Set outaligner_dir parameter depending on active aligner

                if ($verbose) {

                    my $info =
                      q{Set outaligner_dir to: } . $$outaligner_dir_ref;

                    ## Add info to broadcasts
                    push( @$broadcasts_ref, $info );
                }
            }
        }
    }
    if (   ( exists( $aligner{total_active_aligner_count} ) )
        && ( $aligner{total_active_aligner_count} > 1 ) )
    {

        $log->fatal(
            "You have activate more than 1 aligner: "
              . join( ", ", @{ $aligner{active_aligners} } )
              . ". MIP currently only supports 1 aligner per analysis.",
            "\n"
        );
        exit 1;
    }
}

sub collect_read_length {

## Function : Collect read length from an infile
## Returns  : "readLength"
## Arguments: $directory => Directory of file
##          : $file      => File to parse
##          : $read_file => Command used to read file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $file;
    my $read_file_command;

    my $tmpl = {
        directory => {
            defined     => 1,
            required    => 1,
            store       => \$directory,
            strict_type => 1,
        },
        read_file_command => {
            defined     => 1,
            required    => 1,
            store       => \$read_file_command,
            strict_type => 1,
        },
        file =>
          { defined => 1, required => 1, store => \$file, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Prints sequence length and exits
    my $seq_length_regexp =
q?perl -ne 'if ($_!~/@/) {chomp($_);my $seq_length = length($_);print $seq_length;last;}' ?;

    ## Save current direcory
    my $pwd = cwd();

    ## Move to sample_id infile directory
    chdir($directory);

    ## Collect sequence length
    my $ret = `$read_file_command $file | $seq_length_regexp;`;

    ## Move to original directory
    chdir($pwd);

    return $ret;
}

sub check_program_mode {

##check_program_mode

##Function : Check correct value for program mode in MIP.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @allowed_values = ( 0, 1, 2 );

  PROGRAMS:
    foreach my $program ( @{ $parameter_href->{dynamic_parameter}{program} } ) {

        if (
            !(
                any { $_ eq $active_parameter_href->{$program} }
                @allowed_values
            )
          )
        {    #If element is not part of array

            $log->fatal( q{'}
                  . $active_parameter_href->{$program}
                  . q{' Is not an allowed mode for program '--}
                  . $program
                  . q{'. Set to: }
                  . join( "|", @allowed_values ) );
            exit 1;
        }
    }
}

sub check_sample_id_in_parameter_path {

##check_sample_id_in_parameter_path

##Function : Check sample_id provided in hash path parameter is included in the analysis and only represented once
##Returns  : ""
##Tags     : check, sampleids, hash
##Arguments: $active_parameter_href, $sample_ids_ref, $parameter_name
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $sample_ids_ref        => Array to loop in for parameter {REF}
##         : $parameter_names_ref   => Parameter name list {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_ids_ref;
    my $parameter_names_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref
        },
        parameter_names_ref => {
            required => 1,
            defined  => 1,
            default  => [],
            store    => \$parameter_names_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    foreach my $parameter_name (@$parameter_names_ref)
    {    #Lopp through all hash parameters supplied

        my %seen;    #Hash to test duplicate sample_ids later

        foreach my $key ( keys %{ $active_parameter_href->{$parameter_name} } )
        {

            my @parameter_samples =
              split( ",", $active_parameter_href->{$parameter_name}{$key} );

            foreach my $sample_id (@parameter_samples) {

                $seen{$sample_id}++
                  ;    #Increment instance to check duplicates later

                if ( $seen{$sample_id} > 1 ) {    #Check sample_id are unique

                    $log->fatal(
                        "Sample_id: "
                          . $sample_id
                          . " is not uniqe in '-"
                          . $parameter_name . " '"
                          . $key . "="
                          . join( ",", @parameter_samples ),
                        "\n"
                    );
                    exit 1;
                }
            }
        }
        foreach my $sample_id (@$sample_ids_ref) {

            if ( !( any { $_ eq $sample_id } ( keys %seen ) ) )
            {    #If sample_id is not present in parameter_name hash

                $log->fatal(
                    "Could not detect "
                      . $sample_id
                      . " for '--"
                      . $parameter_name
                      . "'. Provided sample_ids are: "
                      . join( ", ", ( keys %seen ) ),
                    "\n"
                );
                exit 1;
            }
        }
    }
}

sub check_sample_id_in_parameter {

## Function : Check sample_id provided in hash parameter is included in the analysis and only represented once
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Holds all parameters {REF}
##          : $parameter_names_ref   => Parameter name list {REF}
##          : $sample_ids_ref        => Array to loop in for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_names_ref;
    my $parameter_href;
    my $sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        parameter_names_ref => {
            required => 1,
            defined  => 1,
            default  => [],
            store    => \$parameter_names_ref
        },
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Loop through all hash parameters supplied
  PARAMETER:
    foreach my $parameter_name ( @{$parameter_names_ref} ) {

        if ( defined $active_parameter_href->{$parameter_name} ) {

            foreach my $sample_id ( @{$sample_ids_ref} ) {

                ## Check that a value exists
                if (
                    !defined(
                        $active_parameter_href->{$parameter_name}{$sample_id}
                    )
                  )
                {

                    next PARAMETER
                      if ( $parameter_href->{$parameter_name}{mandatory} eq
                        q{no} );

                    $log->fatal(
                        "Could not find value for "
                          . $sample_id
                          . " for parameter '--"
                          . $parameter_name . "'",
                        "\n"
                    );
                    exit 1;
                }

                ## If sample_id is not present in parameter_name hash
                if (
                    !(
                        any { $_ eq $sample_id }
                        ( keys %{ $active_parameter_href->{$parameter_name} } )
                    )
                  )
                {

                    $log->fatal(
                        "Could not detect "
                          . $sample_id
                          . " for parameter '--"
                          . $parameter_name
                          . "'. Provided sample_ids for parameter are: "
                          . join(
                            ", ",
                            (
                                keys
                                  %{ $active_parameter_href->{$parameter_name} }
                            )
                          ),
                        "\n"
                    );
                    exit 1;
                }
            }
        }
    }
}

sub get_matching_values_key {

##get_matching_values_key

##Function : Return the key if the hash value and query match
##Returns  : "key pointing to matched value"
##Arguments: $active_parameter_href, $query_value_ref, $parameter_name
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $query_value_ref       => The value to query in the hash {REF}
##         : $parameter_name        => MIP parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $query_value_ref;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        query_value_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$query_value_ref
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %reversed = reverse %{ $active_parameter_href->{$parameter_name} }
      ;    #Values are now keys and vice versa

    if ( exists $reversed{$$query_value_ref} ) {

        return $reversed{$$query_value_ref};
    }
}

sub check_vcfanno_toml {

##check_vcfanno_toml

##Function : Check that the supplied vcfanno toml frequency file match record 'file=' within toml config file
##Returns  : ""
##Arguments: $vcfanno_file_toml, $vcfanno_file_freq
##         : $vcfanno_file_toml => Toml config file
##         : $vcfanno_file_freq => Frequency file recorded inside toml file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vcfanno_file_toml;
    my $vcfanno_file_freq;

    my $tmpl = {
        vcfanno_file_toml => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$vcfanno_file_toml
        },
        vcfanno_file_freq => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$vcfanno_file_freq
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    open( $FILEHANDLE, "<", $vcfanno_file_toml )
      or
      $log->logdie( "Can't open '" . $vcfanno_file_toml . "': " . $! . "\n" );

    while (<$FILEHANDLE>) {

        chomp $_;                          #Remove newline

        if ( $_ =~ /^file="(\S+)"/ ) {

            my $file_path_freq = $1;

            if ( $file_path_freq ne $vcfanno_file_freq ) {

                $log->fatal( "The supplied vcfanno_config_file: "
                      . $vcfanno_file_freq
                      . " does not match record 'file="
                      . $file_path_freq
                      . "' in the sv_vcfanno_config file: "
                      . $vcfanno_file_toml );
                exit 1;
            }
            last;
        }
    }
    close $FILEHANDLE;
}

sub check_snpsift_keys {

##check_snpsift_keys

##Function : Check that the supplied
##Returns  : ""
##Arguments: $snpsift_annotation_files_href, $snpsift_annotation_outinfo_key_href
##         : $snpsift_annotation_files_href       => Snpsift annotation files {REF}
##         : $snpsift_annotation_outinfo_key_href => File and outinfo key to add to vcf {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $snpsift_annotation_files_href;
    my $snpsift_annotation_outinfo_key_href;

    my $tmpl = {
        snpsift_annotation_files_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$snpsift_annotation_files_href
        },
        snpsift_annotation_outinfo_key_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$snpsift_annotation_outinfo_key_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    foreach my $file ( keys %$snpsift_annotation_outinfo_key_href ) {

        unless ( exists( $snpsift_annotation_files_href->{$file} ) ) {

            $log->fatal( "The supplied snpsift_annotation_outinfo_key file: "
                  . $file
                  . " does not match any file in '--snpsift_annotation_files'"
            );
            $log->fatal( "Supplied snpsift_annotation_files files:\n"
                  . join( "\n", keys %$snpsift_annotation_files_href ) );
            exit 1;
        }
    }
}

sub check_key_exists_in_hash {

##Function : Test if key from query hash exists truth hash
##Returns  :
##Arguments: $parameter_name => Parameter name
##         : $truth_href     => Truth hash
##         : $query_href     => Query hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_name;
    my $truth_href;
    my $query_href;

    my $tmpl = {
        parameter_name =>
          { required => 1, defined => 1, store => \$parameter_name },
        truth_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$truth_href
        },
        query_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$query_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

  QUERY_KEY:
    foreach my $key ( keys %{$query_href} ) {

        if ( not exists( $truth_href->{$key} ) ) {

            $log->fatal( $parameter_name
                  . q{ key '}
                  . $key
                  . q{' - Does not exist as module program parameter in MIP} );
            exit 1;
        }
    }
    return;
}

sub check_element_exists_in_hash {

##check_element_exists_in_hash

##Function : Test if element from query array exists truth hash
##Returns  : ""
##Arguments: $truth_href, $queryies, $parameter_name
##         : $truth_href     => Truth hash
##         : $queryies       => Query array
##         : $parameter_name => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $truth_href;
    my $queryies;
    my $parameter_name;

    my $tmpl = {
        truth_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$truth_href
        },
        queryies => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$queryies
        },
        parameter_name =>
          { required => 1, defined => 1, store => \$parameter_name },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    foreach my $element ( @{$queryies} ) {

        if ( !exists( $truth_href->{$element} ) ) {

            $log->fatal( $parameter_name
                  . " element '"
                  . $element
                  . "' - Does not exist as module program parameter in MIP" );
            exit 1;
        }
    }
}

##Investigate potential autodie error
if ( $@ and $@->isa("autodie::exception") ) {

    if ( $@->matches("default") ) {

        say "Not an autodie error at all";
    }
    if ( $@->matches("open") ) {

        say "Error from open";
    }
    if ( $@->matches(":io") ) {

        say "Non-open, IO error.\n";
    }
}
elsif ($@) {

    say "A non-autodie exception.";
}

1;
