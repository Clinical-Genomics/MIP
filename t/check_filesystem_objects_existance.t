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
use Params::Check qw{ check allow last_error };

use FindBin qw{ $Bin };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir };
use Getopt::Long;
use Test::More;
use File::Temp qw{ tempdir tempfile };

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};
Readonly my $COMMA   => q{,};

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
    my %perl_module;

    $perl_module{q{MIP::Script::Utils}} = [qw{ help }];

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Check::Path});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Path qw{ check_filesystem_objects_existance };

diag(   q{Test check_filesystem_objects_existance from Path.pm v}
      . $MIP::Check::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create a temp dir
my $dir = tempdir( CLEANUP => 1 );

## Create a temp file using newly created temp dir
my ( $fh, $file_name ) = tempfile( DIR => $dir );

my %parameter = (
    dir       => $dir,
    file_name => { build_file => 1 },
);

### TEST
## Dirs
my ($exist) = check_filesystem_objects_existance(
    {
        parameter_name => q{dir},
        object_name    => $dir,
        object_type    => q{directory},
    }
);

is( $exist, 1, q{Found directory} );

($exist) = check_filesystem_objects_existance(
    {
        parameter_name => q{dir},
        object_name    => q{does_not_exist},
        object_type    => q{directory},
    }
);

is( $exist, 0, q{No directory} );

## Files
($exist) = check_filesystem_objects_existance(
    {
        parameter_name => q{file_name},
        object_name    => $file_name,
        object_type    => q{file},
    }
);

is( $exist, 1, q{Found file} );

($exist) = check_filesystem_objects_existance(
    {
        parameter_name => q{file_name},
        object_name    => q{does_not_exist},
        object_type    => q{file},
    }
);

is( $exist, 0, q{No file} );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

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

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
