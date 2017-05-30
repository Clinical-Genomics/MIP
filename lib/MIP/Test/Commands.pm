package MIP::Test::Commands;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Test::More;

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = '1.00';

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(generate_call);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

sub generate_call {

##generate_call

##Function : Generate call to module method
##Returns  : "@commands"
##Arguments: $parameter_href, $module_function_cref, $function_base_command
##         : $parameter_href           => Parameters to submit to module method
##         : $required_parameters_href => Required parameters
##         : $module_function_cref     => Module method to test
##         : $function_base_command    => Function base command

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $required_parameters_href;
    my $function_base_command;
    my $module_function_cref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        required_parameters_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$required_parameters_href
        },
        function_base_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$function_base_command
        },
        module_function_cref =>
          { required => 1, defined => 1, store => \$module_function_cref },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    foreach my $parameter ( keys %{$parameter_href} ) {

        my $input_value = $parameter_href->{$parameter}{input};

        ## Store commands from module function
        my @commands;

        if ( %{$required_parameters_href} ) {

            ## Collect required keys and values to generate args
            my @keys = keys %{$required_parameters_href};
            my @values =
              map { ( $required_parameters_href->{$_}{input} ) } @keys;

            ## Collect the parameter keys and values to test
            push @keys,   $parameter;
            push @values, $input_value;

            ## Build arguments
            my @args;
            while ( my ( $key_index, $key ) = each @keys ) {

                push @args, $keys[$key_index], $values[$key_index];
            }

            ## Submit arguments to coderef sub
            @commands = $module_function_cref->( {@args} );

        }
        else {

            ## Submit arguments to coderef sub
            @commands =
              $module_function_cref->( { $parameter => $input_value, } );
        }

        my $SPACE = q{ };

        ## Base command and parameter in question
        my $parameter_to_test = join $SPACE, @commands;

        ## Expected return value from sub call
        my $expected_return = $parameter_href->{$parameter}{expected_output};

        ## Perform test
        ok( $parameter_to_test =~ /$expected_return/,
            'Test argument: ' . $parameter );
    }
    return;
}

1;
