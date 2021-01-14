package MIP::Language::Java;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $EQUALS $SPACE };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ java_core };

}

sub java_core {

## Function : Perl wrapper for writing java recipe to $filehandle. Based on java openjdk version "1.8.0_92".
## Returns  : @commands
## Arguments: $filehandle                => Filehandle to write to
##          : $java_jar                  => The JAR {Optional}
##          : $java_use_large_pages      => Use java large pages {Optional}
##          : $memory_allocation         => Memory allocation for java
##          : $picard_use_barclay_parser => Use legacy CLI parser for picard
##          : $temp_directory            => Redirect tmp files to java temp {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $java_jar;
    my $java_use_large_pages;
    my $memory_allocation;
    my $picard_use_barclay_parser;
    my $temp_directory;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        java_jar => {
            store       => \$java_jar,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        picard_use_barclay_parser => {
            store       => \$picard_use_barclay_parser,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Java core
    # Stores commands depending on input parameters
    my @commands = qw{ java };

    if ($picard_use_barclay_parser) {

        push @commands, q{-Dpicard.useLegacyParser} . $EQUALS . q{false};
    }
    if ($memory_allocation) {

        push @commands, q{-} . $memory_allocation;
    }

    # UseLargePages for requiring large memory pages (cross-platform flag)
    if ($java_use_large_pages) {

        push @commands, q{-XX:-UseLargePages};
    }

    # Temporary directory
    if ($temp_directory) {

        push @commands, q{-Djava.io.tmpdir=} . $temp_directory;
    }

    if ($java_jar) {

        push @commands, q{-jar} . $SPACE . $java_jar;
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
