package MIP::Get::Parameter;

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
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ get_bin_file_path get_capture_kit get_conda_path get_dynamic_conda_path get_module_parameters get_program_parameters get_program_version get_user_supplied_info };
}

## Constants
Readonly my $MINUS_FOUR => -4;
Readonly my $MINUS_ONE  => -1;
Readonly my $MINUS_TWO  => -2;
Readonly my $SPACE      => q{ };

sub get_bin_file_path {

## Function : Get the absolute path to the binary file
## Returns  : $bin_file_path
## Arguments: $active_parameter_href => Hash with active parameters {REF}
##          : $bin_file              => Name of binary file
##          : $conda_path            => Path to conda directory
##          : $environment_href      => Hash with programs and their environments {REF}
##          : $environment_key       => Key to the environment_href [pprogram]

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

        $environment = @{ $environment_href->{$environment_key} }[$MINUS_ONE];
        $bin_file_path =
          catfile( $conda_path, q{envs}, $environment, q{bin}, $bin_file );
    }
    ## Check if main environment in use
    elsif ( $active_parameter_href->{source_main_environment_commands} ) {

        $environment =
          @{ $active_parameter_href->{source_main_environment_commands} }
          [$MINUS_ONE];
        $bin_file_path =
          catfile( $conda_path, q{envs}, $environment, q{bin}, $bin_file );
    }
    ## Assume installed in conda base environment
    else {

        $environment = q{base};
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
        capture_kit => { store => \$capture_kit, strict_type => 1, },
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

    use IPC::Cmd qw{ can_run };

    ## Find path to conda bin
    my $conda_path = can_run(q{conda});

    ## Split ditrs to array
    my @conda_path_dirs = File::Spec->splitdir($conda_path);

    ## Running from conda_environment
    if ( $conda_path_dirs[$MINUS_FOUR] eq q{envs} ) {

        splice @conda_path_dirs, $MINUS_FOUR;
    }
    ## Running from conda base environment
    else {

        splice @conda_path_dirs, $MINUS_TWO;
    }

    ## Return path to conda folder
    return catdir(@conda_path_dirs);
}

sub get_dynamic_conda_path {

## Function : Attempts to find path to directory with binary in conda env
## Returns  : Path to directory
## Arguments: $active_parameters_href => Active parameter hash {REF}
##          : $bin_file               => Bin file to test
##          : $environment_key        => Key to conda environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $bin_file;
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
        environment_key => {
            store => \$environment_key,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Establish path to conda
    if ( not $active_parameter_href->{conda_path} ) {

        $active_parameter_href->{conda_path} = get_conda_path();
    }
    if ( not -d $active_parameter_href->{conda_path} ) {
        return q{Failed to find default conda path};
    }
    my $conda_path = $active_parameter_href->{conda_path};

    ## Get module and program environments in use
    my %environment;
    if ( $active_parameter_href->{program_source_environment_command} ) {
        ## Build hash with "pprogram_name" as keys and "source env command" as value
        @environment{
            keys
              %{ $active_parameter_href->{program_source_environment_command} }
          } =
          values
          %{ $active_parameter_href->{program_source_environment_command} };
    }
    if ( $active_parameter_href->{module_source_environment_command} ) {
        ## Add to environment hash with "pprogram_name" as keys and "source env command" as value
        @environment{
            keys %{ $active_parameter_href->{module_source_environment_command}
            }
          } =
          values %{ $active_parameter_href->{module_source_environment_command}
          };
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
    if ( not $bin_file_path ) {
        return
            q{Failed to find default path for}
          . $SPACE
          . $bin_file
          . $SPACE
          . q{in conda environment}
          . $SPACE
          . $environment;
    }
    if ( not -f $bin_file_path ) {
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

sub get_module_parameters {

## Function : Get core number, time and source environment command
## Returns  : $core_number, $time, @source_environment_cmds
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $mip_program_name      => MIP program name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $mip_program_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        mip_program_name => {
            defined     => 1,
            required    => 1,
            store       => \$mip_program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Initilize variable
    my @source_environment_cmds;

    if (
        exists $active_parameter_href->{module_source_environment_command}
        {$mip_program_name} )
    {

        @source_environment_cmds =
          @{ $active_parameter_href->{module_source_environment_command}
              {$mip_program_name} };
    }
    elsif ( $active_parameter_href->{source_main_environment_commands}
        && @{ $active_parameter_href->{source_main_environment_commands} } )
    {

        @source_environment_cmds =
          @{ $active_parameter_href->{source_main_environment_commands} };
    }
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

    return $core_number, $time, @source_environment_cmds;
}

sub get_program_parameters {

##Function : Get specific source environment command for program
##Returns  : @source_environment_cmds
##Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##         : $mip_program_name      => MIP program name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $mip_program_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        mip_program_name => {
            defined     => 1,
            required    => 1,
            store       => \$mip_program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Initilize variable
    my @source_environment_cmds;

    if (
        exists $active_parameter_href->{program_source_environment_command}
        {$mip_program_name} )
    {

        @source_environment_cmds =
          @{ $active_parameter_href->{program_source_environment_command}
              {$mip_program_name} };
    }
    return @source_environment_cmds;
}

sub get_program_version {

## Function : Get program version by 1. regexp or 2. cmd
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $cmd                   => Command line call
##          : $parameter_name        => Parameter to add version from
##          : $regexp                => Regexp to use for getting version
##          : $sample_info_href      => Info on samples and family hash {REF}

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

        ## Dry run mode for all program
        return if ( $active_parameter_href->{dry_run_all} );

        ## Dry run mode for this program
        return if ( $active_parameter_href->{$parameter_name} eq q{2} );

        my ($version) =
          $active_parameter_href->{$parameter_name} =~ /$regexp/xsm;

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
        sample_origin         => 0,
        supported_capture_kit => 0,
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
