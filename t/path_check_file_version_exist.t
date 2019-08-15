#!/usr/bin/env perl

use Modern::Perl qw{ 2018 };
use warnings qw{ FATAL utf8 };
use autodie;
use 5.026;
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error};
use Cwd;

# Find directory of script
use FindBin qw{ $Bin };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catfile catdir devnull };
use File::Path qw{ make_path remove_tree };
use Getopt::Long;
use Test::More;
use Readonly;
use IPC::Cmd qw{ can_run run };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 0;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA   => q{,};
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

###User Options
GetOptions(
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },    #Display help text
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
        exit;
    },    #Display version number
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
    my %perl_module;

    $perl_module{q{MIP::Script::Utils}} = [qw{ help }];

  PERL_MODULES:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

    ## Modules
    my @modules = (q{MIP::Check::Path});

  MODULES:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Path qw{ check_file_version_exist };

diag(   q{Test check_file_version_exist from Path.pm v}
      . $MIP::Check::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# File path prefix
my $file_prefix = q{test_check_file_version_exist} . $DOT;

# File path suffix
my $file_suffix = $DOT . q{sh};

# File counter
my $file_counter = 0;

# Temporary directory
my $temp_dir = catdir( cwd(), qw{ test_dir .test_check_file_version_exist } );

# Create path
make_path($temp_dir);

my $test_file_path = catfile( $temp_dir, $file_prefix . $file_counter . $file_suffix );

# Open filehandle
open $FILEHANDLE, q{>}, $test_file_path
  or croak( q{Cannot write to '} . $test_file_path . q{' :} . $OS_ERROR . $NEWLINE );

# Create file
print {$FILEHANDLE} q{Test};

close $FILEHANDLE;

## Testing write to file
ok( -e $test_file_path, q{Create test file} );

my ( $file_name, $file_name_tracker ) = check_file_version_exist(
    {
        file_path_prefix => catfile( $temp_dir, $file_prefix ),
        file_path_suffix => $file_suffix,
    }
);

## Test
is( $file_name_tracker, 1, q{File version} );

my $expected_file_counter = 1;

my $expected_filename =
  catfile( $temp_dir, $file_prefix . $expected_file_counter . $file_suffix );
is( $file_name, $expected_filename, q{File name} );

# Clean-up after test
remove_tree( dirname($temp_dir) );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $program_name
##         : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
