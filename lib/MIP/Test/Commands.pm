package MIP::Test::Commands;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw  { catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use List::MoreUtils qw{ zip };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Writefile qw{ test_write_to_file };

## Constants
Readonly my $ERROR_MSG_INDENT => 3;
Readonly my $SPACE => q{ };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ test_function };
}

sub test_function {

## Function : Test module function by generating arguments and testing output
## Returns  : "@commands"
## Arguments: $argument_href              => Parameters to submit to module method
##          : $base_command_index         => index for slicing the $function_base_command_ref
##          : $do_test_base_command       => Perform test of base command
##          : $function_base_commands_ref => Function base commands {REF}
##          : $module_function_cref       => Module method to test
##          : $required_argument_href     => Required arguments
##          : $is_self_testing            => Testing self

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $argument_href;
    my $function_base_commands_ref;
    my $module_function_cref;
    my $required_argument_href;

    ## Default(s)
    my $base_commands_index;
    my $do_test_base_command;
    my $is_self_testing;

    my $tmpl = {
        argument_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$argument_href,
            strict_type => 1,
        },
        base_commands_index => {
            allow       => qr/ ^\d+$ /sxm,
            default     => scalar @{ $arg_href->{function_base_commands_ref} } - 1,
            store       => \$base_commands_index,
            strict_type => 1,
        },
        do_test_base_command => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$do_test_base_command,
            strict_type => 1,
        },
        function_base_commands_ref => {
            defined     => 1,
            default     => [],
            required    => 1,
            store       => \$function_base_commands_ref,
            strict_type => 1,
        },
        module_function_cref => {
            defined  => 1,
            required => 1,
            store    => \$module_function_cref,
        },
        required_argument_href => {
            default     => {},
            store       => \$required_argument_href,
            strict_type => 1,
        },
        is_self_testing => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$is_self_testing,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q(Could not parse arguments!);

  ARGUMENT:
    foreach my $argument ( keys %{$argument_href} ) {

        ### Parameter to test in this loop check if array ref or scalar
        my $input_value;
        my $input_values_ref;

        ## SCALAR
        if ( exists $argument_href->{$argument}{input} ) {

            $input_value = $argument_href->{$argument}{input};
        }
        ## ARRAY
        elsif ( exists $argument_href->{$argument}{inputs_ref} ) {

            $input_values_ref = $argument_href->{$argument}{inputs_ref};
        }

        ## Store commands from module function
        my @commands;

        ## Some functions have mandatory arguments
        if ( %{$required_argument_href} ) {

            my @args;
            if ($input_values_ref) {

                @args = _build_call(
                    {
                        argument               => $argument,
                        input_values_ref       => $input_values_ref,
                        required_argument_href => $required_argument_href,
                    }
                );
            }
            else {

                @args = _build_call(
                    {
                        argument               => $argument,
                        input_value            => $input_value,
                        required_argument_href => $required_argument_href,
                    }
                );
            }

            ## Special case for test of FILEHANDLE. Does not return @commands
            if ( $argument eq q{FILEHANDLE} ) {

                test_write_to_file(
                    {
                        args_ref             => \@args,
                        base_commands_ref    => $function_base_commands_ref,
                        module_function_cref => $module_function_cref,
                    }
                );
            }
            else {

                ## Submit arguments to coderef sub
                @commands = $module_function_cref->( {@args} );
            }
        }
        else {

            ## Special case for test of FILEHANDLE. Does not return @commands
            if ( $argument eq q{FILEHANDLE} ) {

                test_write_to_file(
                    {
                        args_ref             => [ $argument, $input_value ],
                        module_function_cref => $module_function_cref,
                        base_commands_ref    => $function_base_commands_ref,
                    }
                );
            }
            else {

                ## Array
                if ($input_values_ref) {

                    @commands =
                      $module_function_cref->( { $argument => $input_values_ref, } );
                }
                else {

                    ## Submit arguments to coderef sub
                    @commands = $module_function_cref->( { $argument => $input_value, } );

                }
            }
        }

        ## Expected return value from sub call
        my $expected_return = $argument_href->{$argument}{expected_output};

        ### Perform tests

        if (@commands) {

            ## Test function_base_command
            _test_base_command(
                {
                    base_commands_ref    => [ @commands[ 0 .. $base_commands_index ] ],
                    do_test_base_command => $do_test_base_command,
                    expected_base_commands_ref => $function_base_commands_ref,
                    is_self_testing            => $is_self_testing,
                }
            );

            # Remap to hash
            my %is_argument = map { $_ => $_ } @commands;

            if ( exists $is_argument{$expected_return} ) {

                is( $is_argument{$expected_return},
                    $expected_return, q{Argument: } . $argument );
            }
            else {
              TODO: {
                    local $TODO = q{Self testing should fail in test_function.t}, 1,
                      if ($is_self_testing);

                    is(
                        join( $SPACE, @commands ),
                        $expected_return,
                        q{Argument: } . $argument
                    );
                    say {*STDERR} q{#}
                      . $SPACE x $ERROR_MSG_INDENT
                      . q{Command line does not contain expected argument.};
                }
            }
        }
    }
    return;
}

sub _build_call {

## Function : Build arguments to function
## Returns  : "@arguments"
## Arguments: $argument               => Argument key to test
##          : $input_value            => Argument value to test
##          : $input_values_ref       => Argument values to test
##          : $required_argument_href => Required arguments

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $argument;
    my $input_value;
    my $input_values_ref;
    my $required_argument_href;

    my $tmpl = {
        argument => {
            store       => \$argument,
            strict_type => 1,
        },
        input_value => {
            store       => \$input_value,
            strict_type => 1,
        },
        input_values_ref => {
            default     => [],
            store       => \$input_values_ref,
            strict_type => 1,
        },
        required_argument_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$required_argument_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q(Could not parse arguments!);

    ## Collect required keys and values to generate args
    my @keys;
    my @possible_input_names = qw{ input inputs_ref };
    my @values;

  REQUIRED_ARGUMENT:
    foreach my $required_argument ( keys %{$required_argument_href} ) {

        # Add required_argument
        push @keys, $required_argument;

        # Add value
      POSSIBLE_INPUT_NAMES:
        foreach my $input_name (@possible_input_names) {

            ## SCALAR or ARRAY_ref
            if ( exists $required_argument_href->{$required_argument}{$input_name} ) {

                push @values, $required_argument_href->{$required_argument}{$input_name};
            }
        }
    }

    ### Combine the specific and required argument keys and values to test
    ## SCALAR
    if ( $argument && $input_value ) {
        push @keys,   $argument;
        push @values, $input_value;
    }
    ## ARRAY
    elsif ( $argument && @{$input_values_ref} ) {
        push @keys,   $argument;
        push @values, $input_values_ref;
    }

    ## Interleave arrays to build arguments for submission to function
    my @args = zip( @keys, @values );

    return @args;
}

sub _test_base_command {

## Function : Test the function base command. Executable, ".jar" etc.
## Returns  : ""
## Arguments: $base_commands_ref          => Base command(s) {REF}
##          : $do_test_base_command       => Perform test of base command
##          : $expected_base_commands_ref => Expected base command(s) {REF}
##          : $is_self_testing            => Testing self

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $base_commands_ref;
    my $expected_base_commands_ref;

    ## Default(s)
    my $do_test_base_command;
    my $is_self_testing;

    my $tmpl = {
        base_commands_ref => {
            defined     => 1,
            default     => [],
            required    => 1,
            store       => \$base_commands_ref,
            strict_type => 1,
        },
        do_test_base_command => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$do_test_base_command,
            strict_type => 1,
        },
        expected_base_commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$expected_base_commands_ref,
            strict_type => 1,
        },
        is_self_testing => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$is_self_testing,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q(Could not parse arguments!);

    return if ( not $do_test_base_command );

    ## Compare base argument string to the expected one
    if ( ( join $SPACE, @{$base_commands_ref} ) ne
        ( join $SPACE, @{$expected_base_commands_ref} ) )
    {

      TODO: {
            local $TODO = q{Self testing should fail in test_function.t}, 1,
              if ($is_self_testing);

            ## Display base argument difference and exit
            is_deeply(
                $base_commands_ref,         $expected_base_commands_ref,
                'Argument: ' . join $SPACE, @{$expected_base_commands_ref}
            );
        }
    }
    return;
}

1;
