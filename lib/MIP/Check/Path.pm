package MIP::Check::Path;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use autodie;
use List::MoreUtils qw{ any };
use Readonly;

## MIPs lib
use MIP::Constants qw{ $AMPERSAND $CLOSE_BRACKET $OPEN_BRACKET $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.14;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_file_version_exist
      check_future_filesystem_for_directory
    };
}

sub check_file_version_exist {

## Function : Check if a file with with a filename consisting of $file_path_prefix.$file_counter.$file_path_suffix exist. If so bumps the version number and return new file path and version number.
## Returns  : $file_path, $file_name_version
## Arguments: $file_path_prefix => File path prefix
##          : $file_path_suffix => File path suffix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path_prefix;
    my $file_path_suffix;

    my $tmpl = {
        file_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$file_path_prefix,
            strict_type => 1,
        },
        file_path_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$file_path_suffix,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Nr of sbatch scripts with identical filenames i.e. version number
    my $file_name_version = 0;

    my $file_path = $file_path_prefix . $file_name_version . $file_path_suffix;

  FILE_PATHS:
    while ( -e $file_path ) {

        $file_name_version++;

        ## New file_path to test for existence
        $file_path = $file_path_prefix . $file_name_version . $file_path_suffix;
    }
    return ( $file_path, $file_name_version );
}

sub check_future_filesystem_for_directory {

## Function : Build bash script to check if a directory exists and otherwise create it
## Returns  :
## Arguments: $directory_path => Path to check / create

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory_path;

    my $tmpl = {
        directory_path => {
            defined     => 1,
            required    => 1,
            store       => \$directory_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Coreutils qw{gnu_mkdir};

    my $dir_check_command =
        $OPEN_BRACKET
      . $SPACE . q{! -d}
      . $SPACE
      . $directory_path
      . $SPACE
      . $CLOSE_BRACKET
      . $SPACE
      . $AMPERSAND
      . $AMPERSAND
      . $SPACE;
    my @mkdir_commands = gnu_mkdir(
        {
            indirectory_path => $directory_path,
            parents          => 1,
        }
    );

    $dir_check_command .= join $SPACE, @mkdir_commands;

    return $dir_check_command;
}

1;
