package MIP::Main::Download;

use 5.026;
use Carp;
use charnames qw( :full :short );
use Cwd;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use FindBin qw{ $Bin };
use Getopt::Long;
use IO::Handle;
use warnings qw{ FATAL utf8 };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use Time::Piece;
use utf8;

## CPANM
use autodie qw{ open close :all };
use List::Util qw{ any };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use MIP::Check::Parameter
  qw{ check_cmd_config_vs_definition_file check_email_address check_recipe_exists_in_hash check_recipe_mode };
use MIP::Check::Path qw{ check_parameter_files };
use MIP::Cluster qw{ check_max_core_number };
use MIP::Constants
  qw{ $COLON $COMMA $DOT $MIP_VERSION $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Gnu::Coreutils qw{ gnu_mkdir };
use MIP::Language::Shell qw{ create_bash_file };
use MIP::Log::MIP_log4perl qw{ initiate_logger set_default_log4perl_file };
use MIP::Set::Parameter
  qw{ set_config_to_active_parameters set_custom_default_to_active_parameter set_default_to_active_parameter set_cache };
use MIP::Parse::Parameter qw{ parse_download_reference_parameter };
use MIP::Recipes::Pipeline::Download_rd_dna qw{ pipeline_download_rd_dna };
use MIP::Recipes::Pipeline::Download_rd_rna qw{ pipeline_download_rd_rna };
use MIP::Script::Utils qw{ create_temp_dir };
use MIP::Update::Path qw{ update_to_absolute_path };
use MIP::Update::Recipes qw{ update_recipe_mode_with_dry_run_all };

BEGIN {
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables that can be optionally exported
    our @EXPORT_OK = qw{ mip_download };

}

sub mip_download {

## Function : Main script for generating MIP download scripts
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
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

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments};

    ## Transfer to lexical variables
    # Parameters to include in each download run
    my %active_parameter = %{$active_parameter_href};

    # All parameters MIP download knows
    my %parameter = %{$parameter_href};

    ## Get local time
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;
    my $date            = $date_time->ymd;

    # Catches script name and removes ending
    my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{pl} ) );

    # Add specific MIP process
    $script .= $UNDERSCORE . q{download};

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
    }

    ## Set the default Log4perl file using supplied dynamic parameters.
    $active_parameter{log_file} = set_default_log4perl_file(
        {
            cmd_input       => $active_parameter{log_file},
            date            => $date,
            date_time_stamp => $date_time_stamp,
            script          => $script,
        }
    );

    ## Initiate logger
    my $log = initiate_logger(
        {
            file_path => $active_parameter{log_file},
            log_name  => uc q{mip_download},
        }
    );
    $log->info( q{MIP Version: } . $MIP_VERSION );
    $log->info(
        q{Writing log messages to} . $COLON . $SPACE . $active_parameter{log_file} );

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
        my @custom_default_parameters = qw{ reference_dir };

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
                log            => $log,
                parameter_href => \%parameter,
                parameter_name => $parameter_name,
            }
        );
    }

    ## Make sure that we have lower case from user input
    @{ $active_parameter{reference_genome_versions} } =
      map { lc } @{ $active_parameter{reference_genome_versions} };

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

    ## Check email adress syntax and mail host
    check_email_address(
        {
            email => $active_parameter{email},
            log   => $log,
        }
    );

    ## Parameters that have keys as MIP recipe names
    my @parameter_keys_to_check = (qw{ recipe_time recipe_core_number });
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

    ## Check that the module core number do not exceed the maximum per node
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

    ## Adds dynamic aggregate information from definitions to parameter hash
    set_cache(
        {
            aggregates_ref => [
                ## Collects all recipes that MIP can handle
                q{type:recipe},
            ],
            parameter_href => \%parameter,
        }
    );

    ## Check correct value for recipe mode in MIP
    check_recipe_mode(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            parameter_href        => \%parameter,
        }
    );

    ## Update recipe mode depending on dry_run_all flag
    update_recipe_mode_with_dry_run_all(
        {
            active_parameter_href => \%active_parameter,
            dry_run_all           => $active_parameter{dry_run_all},
            recipes_ref           => \@{ $parameter{cache}{recipe} },
        }
    );

    ## Remodel depending on if "--reference" was used or not as the user info is stored as a scalar per reference_id while yaml is stored as arrays per reference_id
    parse_download_reference_parameter(
        { reference_href => \%{ $active_parameter{reference} }, } );

    check_user_reference(
        {
            user_supplied_reference_ref => \%{ $active_parameter{reference} },
            reference_genome_versions_ref =>
              \@{ $active_parameter{reference_genome_versions} },
            reference_ref => \%{ $active_parameter{reference_feature} },
        }
    );

    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    # Download temp directory
    my $temp_dir = create_temp_dir( {} );

    if ( $active_parameter{sbatch_mode} ) {

        $log->info(
q{Will write sbatch install instructions to for sbatch enabled references to individual sbatch scripts}
        );

    }

    # Downloads instruction file
    my $bash_file_path = catfile( cwd(), q{download_reference} . $DOT . q{sh} );

    open $FILEHANDLE, '>', $bash_file_path
      or
      $log->logdie( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . "\n" );

## Create bash file for writing install instructions
    create_bash_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            file_name   => $bash_file_path,
            log         => $log,
            remove_dir  => $temp_dir,
            set_errexit => 1,
            set_nounset => 1,
        }
    );

    say {$FILEHANDLE} q{## Create temp dir};
    gnu_mkdir(
        {
            FILEHANDLE       => $FILEHANDLE,
            indirectory_path => $temp_dir,
            parents          => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    $log->info( q{Will write install instructions to '} . $bash_file_path,
        $SINGLE_QUOTE );

    if ( $active_parameter{sbatch_mode} ) {

        my $pipeline_type = $active_parameter{download_pipeline_type};

        ## Create dispatch table of pipelines
        my %pipeline = (
            rd_dna => \&pipeline_download_rd_dna,
            rd_rna => \&pipeline_download_rd_rna,
        );

        $log->info( q{Pipeline download type: } . $pipeline_type );
        $pipeline{$pipeline_type}->(
            {
                active_parameter_href => \%active_parameter,
                FILEHANDLE            => $FILEHANDLE,
                temp_directory        => $temp_dir,
            }
        );
    }
    else {

        ## Build install references recipe in bash file
        build_reference_install_recipe(
            {
                active_parameter_href => \%active_parameter,
                FILEHANDLE            => $FILEHANDLE,
                temp_directory        => $temp_dir,
            }
        );
    }
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    return;
}

#################
###SubRoutines###
#################

sub build_reference_install_recipe {

## Function : Build install references recipe in bash file
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $FILEHANDLE            => Filehandle to write to
##          : $quiet                 => Be quiet
##          : $temp_directory        => Temporary directory
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;

    ## Default(s)
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
        FILEHANDLE => { defined => 1, required => 1, store => \$FILEHANDLE, },
        quiet      => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbose => {
            default     => $arg_href->{active_parameter_href}{verbose},
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Gnu::Coreutils qw{ gnu_mkdir };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Recipes::Download::Human_reference qw{ download_human_reference };

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger(q{Download_reference});

    my $pwd = cwd();

    ### Download recipes
    ## Create code reference table for download recipes
    my %download_recipe = ( human_reference => \&download_human_reference, );

    say {$FILEHANDLE} q{## Create reference directory};
    gnu_mkdir(
        {
            indirectory_path => $active_parameter_href->{reference_dir},
            parents          => 1,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Since all commands should assume working directory to be the reference directory
    gnu_cd(
        {
            directory_path => $active_parameter_href->{reference_dir},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

  REFERENCE:
    while ( my ( $reference_id, $versions_ref ) =
        each %{ $active_parameter_href->{reference} } )
    {

      REFERENCE_VERSION:
        foreach my $reference_version ( @{$versions_ref} ) {

          GENOME_VERSION:
            foreach my $genome_version (
                @{ $active_parameter_href->{reference_genome_versions} } )
            {

                ## Standardize case
                $genome_version = lc $genome_version;

                my $reference_href =
                  $active_parameter_href->{reference_feature}{$reference_id}
                  {$genome_version}{$reference_version};

                next GENOME_VERSION
                  if (
                    not exists $active_parameter_href->{reference_feature}{$reference_id}
                    {$genome_version} );

                next GENOME_VERSION
                  if (
                    not exists $active_parameter_href->{reference_feature}{$reference_id}
                    {$genome_version}{$reference_version} );

                ## Build file name and path
                my $outfile_name = $reference_href->{outfile};
                my $outfile_path =
                  catfile( $active_parameter_href->{reference_dir}, $outfile_name );

                ## Check if reference already exists in reference directory
                next GENOME_VERSION if ( -f $outfile_path );

                $log->warn( q{Cannot find reference file:} . $outfile_path );
                $log->warn(
                        q{Will try to download: }
                      . $reference_id
                      . q{ version: }
                      . $reference_version,
                );

                get_reference(
                    {
                        FILEHANDLE     => $FILEHANDLE,
                        recipe_name    => $reference_id,
                        reference_dir  => $active_parameter_href->{reference_dir},
                        reference_href => $reference_href,
                        quiet          => $quiet,
                        verbose        => $verbose,
                    }
                );

                ## Writes post processing commands associated with reference e.g. tabix
                write_post_processing_command(
                    {
                        reference_href => $reference_href,
                        FILEHANDLE     => $FILEHANDLE,
                    }
                );
            }
        }
    }

    ## Move back to original
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    return;
}

sub check_user_reference {

## Function : Check that the user supplied reference id and version
## Returns  :
## Arguments: $reference_genome_versions_ref => Reference genome build versions
##          : $reference_ref                 => Defined reference id and version
##          : $user_supplied_reference_ref   => User supplied reference id and version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_genome_versions_ref;
    my $reference_ref;
    my $user_supplied_reference_ref;

    my $tmpl = {
        reference_genome_versions_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_versions_ref,
            strict_type => 1,
        },
        reference_ref => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_ref,
            strict_type => 1,
        },
        user_supplied_reference_ref => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$user_supplied_reference_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Store what has been seen
    my %cache;

  REFERENCE:
    while ( my ( $reference_id, $versions_ref ) = each %{$user_supplied_reference_ref} ) {

        if ( not exists $reference_ref->{$reference_id} ) {

            $log->fatal( q{Cannot find reference key:} . $reference_id );
            exit 1;
        }

      REFERENCE_VERSION:
        foreach my $version ( @{$versions_ref} ) {

          GENOME_VERSION:
            foreach my $reference_genome_version ( @{$reference_genome_versions_ref} ) {

                ## Found match
                if (
                    exists $reference_ref->{$reference_id}{$reference_genome_version}
                    {$version} )
                {

                    $cache{$reference_id}++;
                }
                ## Store mismatch
                push @{ $cache{$version} }, $reference_genome_version;
            }

            ## Require at least one match
            next REFERENCE if ( $cache{$reference_id} );

            $log->fatal(
                q{Cannot find version key: }
                  . $version
                  . q{ reference key:}
                  . $reference_id
                  . q{ genome build version:}
                  . join $COMMA,
                @{ $cache{$version} },
            );
            exit 1;
        }

    }
    return;
}

sub write_post_processing_command {

## Function : Writes post processing commands associated with reference e.g. tabic
## Returns  :
## Arguments: $reference_href => Reference hash {REF}
##          : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_href;
    my $FILEHANDLE;

    my $tmpl = {
        reference_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$reference_href,
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Process reference with commands
    my @processing_commands =
      qw{ outfile_reformat_command outfile_bgzip_command outfile_tabix_command};

  COMMAND:
    foreach my $command (@processing_commands) {

        next COMMAND if ( not exists $reference_href->{$command} );

        ## Command
        say {$FILEHANDLE} $reference_href->{$command}, $NEWLINE;
    }
    return;
}

1;
