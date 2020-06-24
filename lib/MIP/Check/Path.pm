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
    our $VERSION = 1.15;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_future_filesystem_for_directory
    };
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
