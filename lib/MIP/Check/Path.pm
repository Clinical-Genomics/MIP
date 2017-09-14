package MIP::Check::Path;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{ check allow last_error };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_file_version_exist check_dir_path_exist };
}

sub check_file_version_exist {

##check_file_version_exist

##Function : Check if a file with with a filename consisting of $file_path_prefix.$file_counter.$file_path_suffix exist. If so bumps the version number and return new file path and version number.
##Returns  : "$file_path, $file_name_version"
##Arguments: $file_path_prefix, $file_path_suffix
##         : $file_path_prefix => File path prefix
##         : $file_path_suffix => File path suffix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path_prefix;
    my $file_path_suffix;

    my $tmpl = {
        file_path_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path_prefix
        },
        file_path_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path_suffix
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Nr of sbatch scripts with identical filenames i.e. version number
    my $file_name_version = 0;

    my $file_path = $file_path_prefix . $file_name_version . $file_path_suffix;

  FILE_PATHS:
    while ( -e $file_path ) {

        $file_name_version++;

        # New file_path to test for existence
        $file_path = $file_path_prefix . $file_name_version . $file_path_suffix;
    }
    return ( $file_path, $file_name_version );
}

sub check_dir_path_exist {

## check_dir_path_exist

## Function  : Checks if any of the supplied paths exists. Returns an array that holds the existing paths.
## Returns   : @existing_dir_paths
## Arguments : $dir_paths_ref
##           : $dir_paths_ref => Ref to array of supplied paths

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dir_paths_ref;

    my $tmpl = {
        dir_paths_ref => {
            default     => [],
            required    => 1,
            strict_type => 1,
            store       => \$dir_paths_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @existing_dir_paths;

    # Search for directory in supplied paths
    @existing_dir_paths = grep { -d } @{$dir_paths_ref};

    return @existing_dir_paths;
}

1;
