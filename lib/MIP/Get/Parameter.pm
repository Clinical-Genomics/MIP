package MIP::Get::Parameter;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
      qw{ get_capture_kit get_module_parameters get_program_parameters get_user_supplied_info };
}

## Constants
Readonly my $SPACE => q{ };

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

        ## Supported capture kit alias
        if ( defined $supported_capture_kit_href->{default}{$capture_kit} ) {

            return $supported_capture_kit_href->{default}{$capture_kit};
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

        ## Supported capture kit alias
        if ( defined $supported_capture_kit_href->{default}{$capture_kit} ) {

            return $supported_capture_kit_href->{default}{$capture_kit};
        }
        else {
            #Return unchanged capture_kit string

            return $capture_kit;
        }
    }
    return;
}

sub get_module_parameters {

##Function : Get core number, time and source environment command
##Returns  : $core_number, $time, $source_environment_cmd
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
    my $source_environment_cmd;

    if (
        exists $active_parameter_href->{module_source_environment_command}
        {$mip_program_name} )
    {

        $source_environment_cmd =
          $active_parameter_href->{module_source_environment_command}
          {$mip_program_name};
    }
    elsif ( $active_parameter_href->{source_main_environment_commands}
        && @{ $active_parameter_href->{source_main_environment_commands} } )
    {

        $source_environment_cmd = join $SPACE,
          @{ $active_parameter_href->{source_main_environment_commands} };
    }
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

    return $core_number, $time, $source_environment_cmd;
}

sub get_program_parameters {

##Function : Get specific source environment command for program
##Returns  : $source_environment_cmd
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
    my $source_environment_cmd;

    if (
        exists $active_parameter_href->{program_source_environment_command}
        {$mip_program_name} )
    {

        $source_environment_cmd =
          $active_parameter_href->{program_source_environment_command}
          {$mip_program_name};
    }
    return $source_environment_cmd;
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
        analysis_type     => 0,
        exome_target_bed  => 0,
        expected_coverage => 0,
        sample_ids        => 0,
        sample_origin     => 0,
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
