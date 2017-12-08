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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_module_parameters get_program_parameters };
}

## Constants
Readonly my $SPACE => q{ };

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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        mip_program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$mip_program_name
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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        mip_program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$mip_program_name
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

1;
