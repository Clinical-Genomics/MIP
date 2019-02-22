package MIP::Get::Parameter;

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
use List::MoreUtils qw { uniq };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $DOT $EMPTY_STR $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.11;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_bin_file_path
      get_capture_kit
      get_conda_path
      get_dynamic_conda_path
      get_env_method_cmds
      get_gatk_intervals
      get_install_parameter_attribute
      get_package_env_attributes
      get_package_source_env_cmds
      get_pedigree_sample_id_attributes
      get_program_executables
      get_program_version
      get_programs_for_shell_installation
      get_read_group
      get_recipe_parameters
      get_recipe_attributes
      get_user_supplied_info
    };
}

## Constants
Readonly my $MINUS_ONE => -1;
Readonly my $MINUS_TWO => -2;

sub get_bin_file_path {

## Function : Get the absolute path to the binary file
## Returns  : $bin_file_path
## Arguments: $active_parameter_href => Hash with active parameters {REF}
##          : $bin_file              => Name of binary file
##          : $conda_path            => Path to conda directory
##          : $environment_href      => Hash with programs and their environments {REF}
##          : $environment_key       => Key to the environment_href [program]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $bin_file;
    my $conda_path;
    my $environment_href;
    my $environment_key;

    my $tmpl = {
        active_parameter_href => {
            default  => {},
            required => 1,
            store    => \$active_parameter_href,
        },
        bin_file => {
            defined  => 1,
            required => 1,
            store    => \$bin_file,
        },
        conda_path => {
            defined  => 1,
            required => 1,
            store    => \$conda_path,
        },
        environment_href => {
            default  => {},
            required => 1,
            store    => \$environment_href,
        },
        environment_key => {
            store => \$environment_key,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Cwd qw{ abs_path };

    ## Get environment and set test path;
    my $environment;
    my $bin_file_path;

    ## Check environments special case env first
    if ( $environment_key and $environment_href->{$environment_key} ) {

        $environment   = @{ $environment_href->{$environment_key} }[$MINUS_ONE];
        $bin_file_path = catfile( $conda_path, q{envs}, $environment, q{bin}, $bin_file );
    }
    ## Assume installed in conda base environment
    else {

        $environment   = q{base};
        $bin_file_path = catfile( $conda_path, q{bin}, $bin_file );
    }

    ## Return absolute path
    return ( abs_path($bin_file_path), $environment );
}

sub get_capture_kit {

## Function : Return a capture kit depending on user info. If arg->{user_supplied_parameter_switchRef} is set, go a head and add capture kit no matter what the switch was.
## Returns  : "$capture kit", "supported capture kit" or "undef"
## Arguments: $capture_kit                    => Capture kit to add
##          : $supported_capture_kit_href     => The supported capture kits hash {REF}
##          : $user_supplied_parameter_switch => Has user supplied parameter {OPTIONAL}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $capture_kit;
    my $user_supplied_parameter_switch;
    my $supported_capture_kit_href;

    my $tmpl = {
        capture_kit                => { store => \$capture_kit, strict_type => 1, },
        supported_capture_kit_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$supported_capture_kit_href,
            strict_type => 1,
        },
        user_supplied_parameter_switch =>
          { store => \$user_supplied_parameter_switch, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set default or return supplied capture kit
    if ( not defined $user_supplied_parameter_switch ) {

        if ( defined $supported_capture_kit_href->{$capture_kit} ) {

            return $supported_capture_kit_href->{$capture_kit};
        }
        else {
            ## Return unchanged capture_kit string

            return $capture_kit;
        }
    }
    ## Only add if user supplied no info on parameter
    if ( defined $user_supplied_parameter_switch
        and not $user_supplied_parameter_switch )
    {

        if ( defined $supported_capture_kit_href->{$capture_kit} ) {

            return $supported_capture_kit_href->{$capture_kit};
        }
        else {
            #Return unchanged capture_kit string

            return $capture_kit;
        }
    }
    return;
}

sub get_conda_path {

## Function: Get path to conda directory
## Returns : $conda_path
## Arguments: $bin_file               => Bin file to test

    my ($arg_href) = @_;

## Default(s)
    my $bin_file;

    my $tmpl = {
        bin_file => {
            default     => q{conda},
            store       => \$bin_file,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};
    use IPC::Cmd qw{ can_run };

    ## Find path to conda bin
    my $conda_path = can_run($bin_file);

    return if ( not $conda_path );

    ## Split dirs to array
    my @conda_path_dirs = File::Spec->splitdir($conda_path);

    ## Traverse to conda dir from binary
    splice @conda_path_dirs, $MINUS_TWO;

    ## Return path to conda folder
    return catdir(@conda_path_dirs);
}

sub get_dynamic_conda_path {

## Function : Attempts to find path to directory with binary in conda env
## Returns  : Path to directory
## Arguments: $active_parameters_href => Active parameter hash {REF}
##          : $bin_file               => Bin file to test
##          : $conda_bin_file         => Conda bin file name
##          : $environment_key        => Key to conda environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $bin_file;
    my $conda_bin_file;
    my $environment_key;

    my $tmpl = {
        active_parameter_href => {
            default  => {},
            required => 1,
            store    => \$active_parameter_href,
        },
        bin_file => {
            defined  => 1,
            required => 1,
            store    => \$bin_file,
        },
        conda_bin_file => {
            default => q{conda},
            store   => \$conda_bin_file,
        },
        environment_key => {
            store => \$environment_key,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Establish path to conda
    if ( not $active_parameter_href->{conda_path} ) {

        $active_parameter_href->{conda_path} =
          get_conda_path( { bin_file => $conda_bin_file, } );
    }
    if (   not $active_parameter_href->{conda_path}
        or not -d $active_parameter_href->{conda_path} )
    {

        return q{Failed to find default conda path};
    }
    my $conda_path = $active_parameter_href->{conda_path};

    ## Get module and program environments in use
    my %environment;

    ## Load env
    my ( $env_name, $env_method ) = get_package_env_attributes(
        {
            active_parameter_href => $active_parameter_href,
            package_name          => $environment_key,
        }
    );
    if ($env_name) {

        ## Get env load command
        my @env_method_cmds = get_env_method_cmds(
            {
                action     => q{load},
                env_name   => $env_name,
                env_method => $env_method,
            }
        );

        ## Add to environment hash with "recipe_name" as keys and "source env command" as value
        $environment{$environment_key} = [@env_method_cmds];
    }

    ## Get the bin file path
    my ( $bin_file_path, $environment ) = get_bin_file_path(
        {
            active_parameter_href => $active_parameter_href,
            bin_file              => $bin_file,
            conda_path            => $conda_path,
            environment_href      => \%environment,
            environment_key       => $environment_key,
        }
    );

    ## Test if path exists
    if (   not $bin_file_path
        or not -f $bin_file_path )
    {
        return
            q{Failed to find default path for}
          . $SPACE
          . $bin_file
          . $SPACE
          . q{in conda environment}
          . $SPACE
          . $environment;
    }

    ## Get directory path
    my @bin_path_dirs = File::Spec->splitdir($bin_file_path);
    pop @bin_path_dirs;

    return catdir(@bin_path_dirs);
}

sub get_env_method_cmds {

## Function : Get the standard load and unload env command for environment method
## Returns  : @env_method_cmds
## Arguments: $action     => What to do with the environment
##          : $env_method => Method used to load environment
##          : $env_name   => Name of environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $action;
    my $env_method;
    my $env_name;

    my $tmpl = {
        action => {
            allow       => [qw{ load unload }],
            defined     => 1,
            required    => 1,
            store       => \$action,
            strict_type => 1,
        },
        env_method => {
            allow       => [qw{ conda }],
            defined     => 1,
            required    => 1,
            store       => \$env_method,
            strict_type => 1,
        },
        env_name => {
            required    => 1,
            store       => \$env_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };

    my %method_cmd = (
        conda => {
            load   => [ ( conda_activate(   { env_name => $env_name, } ), ) ],
            unload => [ ( conda_deactivate( {} ), ) ],
        },
    );

    return ( @{ $method_cmd{$env_method}{$action} } );
}

sub get_gatk_intervals {

## Function : Generate and return interval hash
## Returns  : %gatk_intervals
## Argumetns: $analysis_type         => Analysis type
##          : $contigs_ref           => Contigs to split in file
##          : $exome_target_bed_href => Exome target bed files lnked to sample ids
##          : $file_ending           => File ending to add {Optional}
##          : $FILEHANDLE            => Filehandle to write to
##          : $log                   => Log
##          : $max_cores_per_node    => Maximum core per node
##          : $outdirectory          => Outdirectory
##          : $reference_dir         => MIP reference directory
##          : $sample_id             => Sample_id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type;
    my $contigs_ref;
    my $exome_target_bed_href;
    my $file_ending;
    my $FILEHANDLE;
    my $log;
    my $outdirectory;
    my $reference_dir;
    my $sample_id;

    ## Default(s)
    my $max_cores_per_node;

    my $tmpl = {
        analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_type,
            strict_type => 1,
        },
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        exome_target_bed_href => {
            default     => {},
            store       => \$exome_target_bed_href,
            strict_type => 1,
        },
        file_ending => {
            store       => \$file_ending,
            strict_type => 1,
        },
        FILEHANDLE => { store => \$FILEHANDLE, },
        log        => {
            store => \$log,
        },
        max_cores_per_node => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        outdirectory => {
            store       => \$outdirectory,
            strict_type => 1,
        },
        reference_dir => {
            store       => \$reference_dir,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_exom_target_bed_file };
    use MIP::File::Interval qw{ generate_contig_interval_file };

    ## Store gatk interval
    my %gatk_intervals;

    if ( $analysis_type eq q{wes} ) {

        my $exome_target_bed_file = get_exom_target_bed_file(
            {
                exome_target_bed_href => $exome_target_bed_href,
                file_ending           => $file_ending,
                log                   => $log,
                sample_id             => $sample_id,
            }
        );

        ## Generate contig specific interval_list and return hash with paths
        %gatk_intervals = generate_contig_interval_file(
            {
                contigs_ref           => $contigs_ref,
                exome_target_bed_file => $exome_target_bed_file,
                FILEHANDLE            => $FILEHANDLE,
                max_cores_per_node    => $max_cores_per_node,
                outdirectory          => $outdirectory,
                reference_dir         => $reference_dir,
            }
        );
    }
    else {
        ## Key-value pairs are identical for WGS/WTS
        %gatk_intervals = map { $_ => [$_] } @{$contigs_ref};
    }

    return %gatk_intervals;
}

sub get_install_parameter_attribute {

## Function : Return parameter attribute from hash
## Returns  :
## Arguments: $parameter_href => Holds all parameters {REF}
##          : $parameter_name => Name of key to return

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Make sure that the supplied key exists
    croak(qq{Could not find parameter_name key: '$parameter_name' in hash})
      if ( not exists $parameter_href->{$parameter_name} );

    ## Hash attribute
    if ( ref $parameter_href->{$parameter_name} eq q{HASH} ) {

        return %{ $parameter_href->{$parameter_name} };
    }
    ## ARRAY attribute
    if ( ref $parameter_href->{$parameter_name} eq q{ARRAY} ) {

        return @{ $parameter_href->{$parameter_name} };
    }
    ## Scalar attribute
    return $parameter_href->{$parameter_name};
}

sub get_package_env_attributes {

## Function : Get environment name and method for package (recipe, program or MIP)
## Returns  : $env_name, $env_method or "undef"
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $package_name          => Package name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $package_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        package_name => {
            defined     => 1,
            required    => 1,
            store       => \$package_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  ENV:
    foreach my $env_name ( keys %{ $active_parameter_href->{load_env} } ) {

        ## Found recipe within env
        if ( exists $active_parameter_href->{load_env}{$env_name}{$package_name} ) {

            ## Unpack
            my $env_method = $active_parameter_href->{load_env}{$env_name}{method};
            return $env_name, $env_method;
        }
    }
    return;
}

sub get_package_source_env_cmds {

## Function : Get package source environment commands
## Returns  : @source_environment_cmds
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $package_name          => Package name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $package_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        package_name => {
            defined     => 1,
            required    => 1,
            store       => \$package_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Initilize variable
    my @source_environment_cmds;

    my ( $env_name, $env_method ) = get_package_env_attributes(
        {
            active_parameter_href => $active_parameter_href,
            package_name          => $package_name,
        }
    );

    ## Could not find recipe within env
    if ( not $env_name ) {

        ## Fall back to MIPs MAIN env
        ( $env_name, $env_method ) = get_package_env_attributes(
            {
                active_parameter_href => $active_parameter_href,
                package_name          => q{mip},
            }
        );
    }
    ## Prior to env command special case
    ## for recipes needing addtional processing
    my $prior_to_load_cmd = $active_parameter_href->{load_env}{$env_name}{$package_name};
    if ($prior_to_load_cmd) {

        push @source_environment_cmds, $prior_to_load_cmd;
    }

    ## Get env load command
    my @env_method_cmds = get_env_method_cmds(
        {
            action     => q{load},
            env_name   => $env_name,
            env_method => $env_method,
        }
    );
    push @source_environment_cmds, @env_method_cmds;

    return @source_environment_cmds;
}

sub get_pedigree_sample_id_attributes {

## Function : Get pedigree sample id attribute
## Returns  : $attribute
## Arguments: $attribute        => Attribute key
##          : $sample_id        => Sample id to get attribute for
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        attribute => {
            allow => [
                qw{ analysis_type capture_kit expected_coverage father mother phenotype sample_id sample_name is_from_sample sex }
            ],
            defined     => 1,
            required    => 1,
            store       => \$attribute,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    ## Make sure that the supplied key exists
    croak(qq{Could not find sample_info_name key: '$attribute' in hash})
      if ( not exists $sample_info_href->{sample}{$sample_id}{$attribute} );

    return $sample_info_href->{sample}{$sample_id}{$attribute};
}

sub get_program_executables {

## Function : Get the parameter file program executables per recipe
## Returns  : uniq @program_executables
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

    ## Get all program executables
    my @program_executables;

    my $err_msg = q{No keys 'cache' and 'recipe' in parameter hash};
    croak($err_msg) if ( not exists $parameter_href->{cache}{recipe} );

  RECIPE:
    foreach my $recipe ( @{ $parameter_href->{cache}{recipe} } ) {

        if ( exists $parameter_href->{$recipe}{program_executables} ) {

            push @program_executables,
              @{ $parameter_href->{$recipe}{program_executables} };
        }
    }
    ## Make unique and return
    return uniq(@program_executables);
}

sub get_programs_for_shell_installation {

## Function  : Get the programs that are to be installed via SHELL
## Returns   : @shell_programs
## Arguments : $conda_programs_href        => Hash with conda progrmas {REF}
##           : $log                        => Log
##           : $prefer_shell               => Path to conda environment
##           : $shell_install_programs_ref => Array with programs selected for shell installation {REF}
##           : $shell_programs_href        => Hash with shell programs {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_programs_href;
    my $log;
    my $prefer_shell;
    my $shell_install_programs_ref;
    my $shell_programs_href;

    my $tmpl = {
        conda_programs_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$conda_programs_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        prefer_shell => {
            allow       => [ undef, 0, 1 ],
            required    => 1,
            store       => \$prefer_shell,
            strict_type => 1,
        },
        shell_install_programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$shell_install_programs_ref,
            strict_type => 1,
        },
        shell_programs_href => {
            default  => {},
            required => 1,
            store    => \$shell_programs_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ intersect array_minus unique };

    return if not keys %{$shell_programs_href};

    my @shell_programs = keys %{$shell_programs_href};
    my @conda_programs = keys %{$conda_programs_href};

    if ($prefer_shell) {

        # Only get the selected programs otherwise leave the array unaltered
        if ( @{$shell_install_programs_ref} ) {

            # Get the intersect between the two arrays
            @shell_programs =
              intersect( @shell_programs, @{$shell_install_programs_ref} );
        }
    }
    elsif ( @{$shell_install_programs_ref} ) {

        # Get elements in @shell_programs that are not part of the conda hash
        my @shell_only_programs = array_minus( @shell_programs, @conda_programs );

        # Add the selected program(s) and remove possible duplicates
        @shell_programs = unique( @shell_only_programs, @{$shell_install_programs_ref} );
    }
    else {
        # If no shell preferences only add programs lacking conda counterpart
        @shell_programs = array_minus( @shell_programs, @conda_programs );
    }

    return @shell_programs;
}

sub get_program_version {

## Function : Get program version by 1. regexp or 2. cmd
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $cmd                   => Command line call
##          : $parameter_name        => Parameter to add version from
##          : $regexp                => Regexp to use for getting version
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $cmd;
    my $parameter_name;
    my $regexp;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        cmd => {
            defined     => 1,
            required    => 1,
            store       => \$cmd,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        regexp => {
            defined     => 1,
            required    => 1,
            store       => \$regexp,
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

    use MIP::Unix::System qw{ system_cmd_call };

    if ( exists $active_parameter_href->{$parameter_name}
        && $active_parameter_href->{$parameter_name} )
    {

        ## Dry run mode for all recipes
        return if ( $active_parameter_href->{dry_run_all} );

        ## Dry run mode for this recipe
        return if ( $active_parameter_href->{$parameter_name} eq q{2} );

        my ($version) = $active_parameter_href->{$parameter_name} =~ /$regexp/xsm;

        # If not set - fall back on actually calling program
        if ( not $version ) {

            my %return = system_cmd_call( { command_string => $cmd, } );
            if ( $return{output}[0] ) {

                chomp( $version = $return{output}[0] );
            }
        }
        return $version;
    }
    return;
}

sub get_read_group {

## Function : Builds hash with read group headers
## Returns  : %read_group
## Arguments: $infile_prefix    => Name of Fastq file minus read direction information
##          : $platform         => Sequencing platform
##          : $sample_id        => Sample ID
##          : $sample_info_href => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_prefix;
    my $platform;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        infile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$infile_prefix,
            strict_type => 1,
        },
        platform => {
            defined     => 1,
            required    => 1,
            store       => \$platform,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    ## Traverse down $sample_info_href
    my %fastq_file =
      %{ $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {read_direction_file}{ $infile_prefix . q{_1} } };

    ## RG hash
    my %rg;

    ## Add ID
    $rg{id} = $infile_prefix;

    ## Add platform unit
    $rg{pu} =
        $fastq_file{flowcell}
      . $DOT
      . $fastq_file{lane}
      . $DOT
      . $fastq_file{sample_barcode};

    ## Add sample
    $rg{sm} = $sample_id;

    ## Add platform
    $rg{pl} = $platform;

    ## Add molecular library (Dummy value since the actual LB isn't available)
    $rg{lb} = $sample_id;

    return %rg;
}

sub get_recipe_attributes {

## Function : Return recipe attributes
## Returns  : $attribute | %attribute
## Arguments: $attribute      => Attribute key
##          : $parameter_href => Holds all parameters
##          : $recipe_name    => Recipe name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $parameter_href;
    my $recipe_name;

    my $tmpl = {
        attribute => {
            store       => \$attribute,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            required    => 1,
            defined     => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get attribute value
    if ( defined $attribute && $attribute ) {

        return $parameter_href->{$recipe_name}{$attribute};
    }

    ## Get recipe attribute hash
    return %{ $parameter_href->{$recipe_name} };
}

sub get_recipe_parameters {

## Function : Get core number, time and source environment command
## Returns  : $core_number, $time, @source_environment_cmds
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $recipe_name           => Recipe name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $recipe_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Initilize variable
    my @source_environment_cmds = get_package_source_env_cmds(
        {
            active_parameter_href => $active_parameter_href,
            package_name          => $recipe_name,
        }
    );

    my $core_number = $active_parameter_href->{recipe_core_number}{$recipe_name};
    my $time        = $active_parameter_href->{recipe_time}{$recipe_name};

    return $core_number, $time, @source_environment_cmds;
}

sub get_user_supplied_info {

## Function : Detect if user supplied info on parameters otherwise collected from pedigree
## Returns  : %user_supply_switch - Hash where 1=user input and 0=no user input
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Define what should be checked
    my %user_supply_switch = (
        analysis_type         => 0,
        exome_target_bed      => 0,
        expected_coverage     => 0,
        sample_ids            => 0,
        is_from_sample        => 0,
        supported_capture_kit => 0,
        time_point            => 0,
    );

    ## Detect user supplied info
  USER_PARAMETER:
    foreach my $parameter ( keys %user_supply_switch ) {

        ## If hash and supplied
        if ( ref $active_parameter_href->{$parameter} eq q{HASH}
            && keys %{ $active_parameter_href->{$parameter} } )
        {

            $user_supply_switch{$parameter} = 1;
        }
        elsif ( ref $active_parameter_href->{$parameter} eq q{ARRAY}
            && @{ $active_parameter_href->{$parameter} } )
        {
            ## If array and supplied
            $user_supply_switch{$parameter} = 1;
        }
        elsif ( defined $active_parameter_href->{$parameter}
            and not ref $active_parameter_href->{$parameter} )
        {

            ## If scalar and supplied
            $user_supply_switch{$parameter} = 1;
        }
        else {

            ## No user defined input for parameter
            $user_supply_switch{$parameter} = 0;
        }
    }
    return %user_supply_switch;
}

1;
