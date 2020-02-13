#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catfile catdir devnull };
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use List::Util qw(any);
use Modern::Perl;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Script::Utils qw(help);
use MIP::Test::Writefile qw(test_write_to_file);

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Language::Shell});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Language::Shell qw(build_shebang);

diag(   q{Test build_shebang from Shell.pm v}
      . $MIP::Language::Shell::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $separator = q{\n};

## Base arguments
my $batch_shebang = q{#!};

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => $batch_shebang,
    },
);

my $bash_bin_path =
  catfile( dirname( dirname( devnull() ) ), qw(usr bin env bash) );

## Specific arguments
my %argument = (
    bash_bin_path => {
        input           => $bash_bin_path,
        expected_output => $batch_shebang . $bash_bin_path . q{ --login},
    },
    invoke_login_shell => {
        input           => 1,
        expected_output => $batch_shebang . $bash_bin_path . q{ --login},
    },
);

my @commands = build_shebang(
    {
        bash_bin_path      => $argument{bash_bin_path}{input},
        invoke_login_shell => $argument{invoke_login_shell}{input},
    }
);

## Testing return of commands
foreach my $key ( keys %argument ) {

    # Alias expected output
    my $expected_output = $argument{$key}{expected_output};

    ok( ( any { $_ eq $expected_output } @commands ), q{Argument: } . $key );
}

## Testing write to file

# Fake arguments
my @args = (
    bash_bin_path => $bash_bin_path,
    filehandle    => undef,
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&build_shebang;

my @function_base_commands = ( $batch_shebang . $bash_bin_path );

test_write_to_file(
    {
        args_ref             => \@args,
        module_function_cref => $module_function_cref,
        base_commands_ref    => \@function_base_commands,
        separator            => $separator,
    }
);

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $recipe_name
##         : $recipe_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $recipe_name;

    my $tmpl = {
        recipe_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$recipe_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    return <<"END_USAGE";
 $recipe_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
