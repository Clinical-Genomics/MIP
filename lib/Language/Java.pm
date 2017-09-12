package Language::Java;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Readonly;

## MIPs lib/
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

## Constants
Readonly my $SPACE => q{ };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.0.1;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ java_core };

}

sub java_core {

## java_core

## Function : Perl wrapper for writing java recipe to $FILEHANDLE. Based on java openjdk version "1.8.0_92".
## Returns  : @commands
## Arguments: $memory_allocation, $FILEHANDLE, $temp_directory, $java_jar, $java_use_large_pages
##         : $memory_allocation    => Memory allocation for java
##         : $FILEHANDLE           => Filehandle to write to
##         : $temp_directory       => Redirect tmp files to java temp {Optional}
##         : $java_jar             => The JAR
##         : $java_use_large_pages => Use java large pages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $FILEHANDLE;
    my $java_use_large_pages;
    my $java_jar;
    my $temp_directory;

    my $tmpl = {
        memory_allocation => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$memory_allocation
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
        temp_directory => {
            strict_type => 1,
            store       => \$temp_directory
        },
        java_jar => {
            strict_type => 1,
            store       => \$java_jar
        },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Java core
    # Stores commands depending on input parameters
    my @commands = q{java};

    push @commands, q{-} . $memory_allocation;

    # UseLargePages for requiring large memory pages (cross-platform flag)
    if ($java_use_large_pages) {
        push @commands, q{-XX:-UseLargePages};
    }

    # Temporary directory
    if ( defined $temp_directory ) {
        push @commands, q{-Djava.io.tmpdir=} . $temp_directory;
    }

    if ( defined $java_jar ) {
        push @commands, q{-jar} . $SPACE . $java_jar;
    }

    if ($FILEHANDLE) {
        unix_write_to_file(
            {
                commands_ref => \@commands,
                separator    => $SPACE,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }

    return @commands;
}

1;
