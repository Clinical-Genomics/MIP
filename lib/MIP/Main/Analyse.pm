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
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use POSIX;
use Time::Piece;
use utf8;
use warnings qw{ FATAL utf8 };

## Third party module(s)
use autodie qw{ open close :all };
use IPC::System::Simple;
use Modern::Perl qw{ 2018 };
use Path::Iterator::Rule;
use Readonly;

## MIPs lib/
use MIP::Active_parameter qw{
  check_load_env_packages
  check_recipe_mode
  get_not_allowed_temp_dirs
  parse_recipe_resources
  set_gender_sample_ids
  set_parameter_reference_dir_path
  update_recipe_mode_with_dry_run_all
  update_to_absolute_path };
use MIP::Analysis qw{ check_analysis_type_to_pipeline get_overall_analysis_type };
use MIP::Config qw{ parse_config };
use MIP::Constants qw{ $DOT $EMPTY_STR $MIP_VERSION $NEWLINE $SINGLE_QUOTE $SPACE $TAB };
use MIP::Contigs qw{ set_contigs };
use MIP::Environment::User qw{ check_email_address };
use MIP::File_info qw{ set_dict_contigs set_human_genome_reference_features };
use MIP::File::Path qw{ check_allowed_temp_directory };
use MIP::Io::Recipes qw{ build_file_prefix_tag };
use MIP::Log::MIP_log4perl qw{ get_log };
use MIP::Parameter qw{
  check_recipe_vs_binary_name
  get_cache
  parse_parameter_files
  parse_reference_path
  set_cache
  set_cache_program_executables
  set_default
};
use MIP::Pedigree qw{ create_fam_file
  get_is_trio
  parse_pedigree
};
use MIP::Pipeline qw{ run_analyse_pipeline };
use MIP::Processmanagement::Processes qw{ write_job_ids_to_file };
use MIP::Recipes::Parse qw{ parse_recipes parse_start_with_recipe };
use MIP::Reference qw{ check_human_genome_file_endings };
use MIP::Sample_info qw{
  reload_previous_pedigree_info
  set_no_dry_run_parameters
  write_sample_info_to_file };
use MIP::Store qw{ store_files };
use MIP::Validate::Case qw{ check_sample_ids };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.60;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_analyse };
}

## Constants
Readonly my %RECIPE_PARAMETERS_TO_CHECK => (
    keys => [
        qw{ recipe_core_number recipe_memory recipe_time
          set_recipe_core_number set_recipe_memory set_recipe_time }
    ],
    elements => [qw{ associated_recipe decompose_normalize_references }],
);
Readonly my $MINUS_ONE => -1;

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

    # Holds all active parameters values for broadcasting
    my @broadcasts;

#### Script parameters

## Add date_time_stamp for later use in log and qc_metrics yaml file
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;
    my $date            = $date_time->ymd;

    # Catches script name and removes ending
    my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{pl} ) );

#### Set program parameters

    my ( %job_id, %sample_info );

#### Staging Area
### Get and/or set input parameters

## Change relative path to absolute path for parameter with "update_path: absolute_path" in config
    update_to_absolute_path(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

### Config
    parse_config(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

## Get log object and set log file in active parameters unless already set from cmd
    my $log = get_log(
        {
            active_parameter_href => $active_parameter_href,
            date                  => $date,
            date_time_stamp       => $date_time_stamp,
            log_name              => uc q{mip_analyse},
            script                => $script,
        }
    );

## Write MIP VERSION and log file path
    $log->info( q{MIP Version: } . $MIP_VERSION );
    $log->info( q{Script parameters and info from are saved in file: }
          . $active_parameter_href->{log_file} );

## Pedigree
    parse_pedigree(
        {
            active_parameter_href => $active_parameter_href,
            pedigree_file_path    => $active_parameter_href->{pedigree_file},
            parameter_href        => $parameter_href,
            sample_info_href      => \%sample_info,
        }
    );

    ## Detect if all samples has the same sequencing type and return consensus if reached
    $parameter_href->{cache}{consensus_analysis_type} = get_overall_analysis_type(
        {
            analysis_type_href => $active_parameter_href->{analysis_type},
        }
    );

## Set default from parameter hash to active_parameter for uninitilized parameters
    set_default(
        {
            active_parameter_href => $active_parameter_href,
            custom_default_parameters_ref =>
              $parameter_href->{custom_default_parameters}{default},
            parameter_href => $parameter_href,
        }
    );

    my $consensus_analysis_type = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{consensus_analysis_type},
        }
    );

## Update path for supplied reference(s) associated with parameter that should
## reside in the mip reference directory to full path
    set_parameter_reference_dir_path(
        {
            active_parameter_href => $active_parameter_href,
            parameter_name        => q{human_genome_reference},
        }
    );

## Detect version and source of the human_genome_reference: Source (hg19 or GRCh) and check compression status
    set_human_genome_reference_features(
        {
            file_info_href => $file_info_href,
            human_genome_reference =>
              basename( $active_parameter_href->{human_genome_reference} ),
            parameter_href => $parameter_href,
        }
    );

## Reference in MIP reference directory
    parse_reference_path(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

### Checks
    check_analysis_type_to_pipeline(
        {
            pipeline      => lc _parent_module( {} ),
            analysis_type => $consensus_analysis_type,
        }
    );

## Parse existence of files and directories
    parse_parameter_files(
        {
            active_parameter_href   => $active_parameter_href,
            consensus_analysis_type => $consensus_analysis_type,
            parameter_href          => $parameter_href,
        }
    );

## Updates sample_info hash with previous run pedigree info
    reload_previous_pedigree_info(
        {
            sample_info_href      => \%sample_info,
            sample_info_file_path => $active_parameter_href->{sample_info_file},
        }
    );

## Special case since dict is created with .fastq removed
## Check the existance of associated human genome files
    check_human_genome_file_endings(
        {
            human_genome_reference_file_endings_ref =>
              $file_info_href->{human_genome_reference_file_endings},
            human_genome_reference_path =>
              $active_parameter_href->{human_genome_reference},
            parameter_href => $parameter_href,
            parameter_name => q{human_genome_reference_file_endings},
        }
    );

## Set sequence contigs used in analysis from human genome sequence dict file
    my $dict_file_path = catfile( $active_parameter_href->{reference_dir},
        $file_info_href->{human_genome_reference_name_prefix} . $DOT . q{dict} );

    set_dict_contigs(
        {
            dict_file_path => $dict_file_path,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
        }
    );

## Detect case constellation based on pedigree file
    $parameter_href->{cache}{trio} = get_is_trio(
        {
            active_parameter_href => $active_parameter_href,
            sample_info_href      => \%sample_info,
        }
    );

## Check email adress syntax and mail host
    check_email_address(
        {
            email => $active_parameter_href->{email},
        }
    );

## Check that the temp directory value is allowed
    my @is_not_allowed_temp_dirs =
      get_not_allowed_temp_dirs( { active_parameter_href => $active_parameter_href, } );
    check_allowed_temp_directory(
        {
            not_allowed_paths_ref => \@is_not_allowed_temp_dirs,
            temp_directory        => $active_parameter_href->{temp_directory},
        }
    );

## Parameters that have keys or elements as MIP recipe names
    parse_recipes(
        {
            active_parameter_href   => $active_parameter_href,
            parameter_href          => $parameter_href,
            parameter_to_check_href => \%RECIPE_PARAMETERS_TO_CHECK,
        }
    );

    ## Check core number requested against environment provisioned
    parse_recipe_resources( { active_parameter_href => $active_parameter_href, } );

    ## Check that the case_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
    check_sample_ids(
        {
            case_id        => $active_parameter_href->{case_id},
            sample_ids_ref => $active_parameter_href->{sample_ids},
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
            parameter_href => $parameter_href,
        }
    );

    set_cache_program_executables( { parameter_href => $parameter_href, } );

    ## Check correct value for recipe mode in MIP
    check_recipe_mode(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

    ## Check that package name name are included in MIP as either "mip", "recipe" or "program"
    check_load_env_packages(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

    ## Check that recipe name and program name are not identical
    check_recipe_vs_binary_name(
        {
            parameter_href   => $parameter_href,
            recipe_names_ref => $parameter_href->{cache}{recipe},
        }
    );

    parse_start_with_recipe(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        },
    );

## Update recipe mode depending on dry_run_all flag
    update_recipe_mode_with_dry_run_all(
        {
            active_parameter_href => $active_parameter_href,
            dry_run_all           => $active_parameter_href->{dry_run_all},
            recipes_ref           => $parameter_href->{cache}{recipe},
        }
    );

    ## Set the gender(s) included in current analysis and count them
    set_gender_sample_ids(
        {
            active_parameter_href => $active_parameter_href,
            sample_info_href      => \%sample_info,
        }
    );

### Contigs
## Set contig prefix and contig names depending on reference used
    set_contigs(
        {
            file_info_href => $file_info_href,
            version        => $file_info_href->{human_genome_reference_version},
        }
    );

## Creates all fileendings as the samples is processed depending on the chain of modules activated
    my @order_recipes = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{order_recipes_ref},
        }
    );
    build_file_prefix_tag(
        {
            active_parameter_href => $active_parameter_href,
            case_id               => $active_parameter_href->{case_id},
            file_info_href        => $file_info_href,
            order_recipes_ref     => \@order_recipes,
            parameter_href        => $parameter_href,
            sample_ids_ref        => $active_parameter_href->{sample_ids},
        }
    );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            case_id          => $active_parameter_href->{case_id},
            execution_mode   => q{system},
            fam_file_path    => $active_parameter_href->{pedigree_fam_file},
            sample_ids_ref   => $active_parameter_href->{sample_ids},
            sample_info_href => \%sample_info,
        }
    );

############
####MAIN####
############

    set_no_dry_run_parameters(
        {
            is_dry_run_all   => $active_parameter_href->{dry_run_all},
            analysis_date    => $date_time_stamp,
            mip_version      => $MIP_VERSION,
            sample_info_href => \%sample_info,
        }
    );

    run_analyse_pipeline(
        {
            active_parameter_href   => $active_parameter_href,
            broadcasts_ref          => \@broadcasts,
            consensus_analysis_type => $consensus_analysis_type,
            file_info_href          => $file_info_href,
            job_id_href             => \%job_id,
            order_parameters_ref    => $order_parameters_ref,
            order_recipes_ref       => $parameter_href->{cache}{order_recipes_ref},
            parameter_href          => $parameter_href,
            sample_info_href        => \%sample_info,
        }
    );

    ## Write QC for recipes used in analysis
    # Write sample info to yaml file
    write_sample_info_to_file(
        {
            sample_info_file => $active_parameter_href->{sample_info_file},
            sample_info_href => \%sample_info,
        }
    );

    write_job_ids_to_file(
        {
            case_id           => $active_parameter_href->{case_id},
            job_id_href       => \%job_id,
            job_ids_file_path => catfile(
                $active_parameter_href->{outdata_dir},
                q{slurm_job_ids} . $DOT . q{yaml}
            ),
        }
    );

    store_files(
        {
            active_parameter_href => $active_parameter_href,
            sample_info_href      => \%sample_info,
        }
    );
    return;
}

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

sub _parent_module {

## Function : Returns the name of the module that called this one
## Returns  : $parent_module
## Arguments:

    ## Get full path to module
    my $parent_module = ( caller 1 )[0];

    ## Isolate module
    $parent_module = ( split /::/xms, $parent_module )[$MINUS_ONE];

    return $parent_module;
}

1;
