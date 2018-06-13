package MIP::Program::Base::Picardtools;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Params::Check qw{ check allow last_error };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ picardtools_base };
}

## Constants
Readonly my $SPACE => q{ };

sub picardtools_base {

## Function : Perl wrapper for picardtools base. Based on Picardtools v2.9.2-SNAPSHOT
## Returns  : @commands
## Arguments: $commands_ref       => List of commands added earlier
##          : $FILEHANDLE         => Filehandle to write to
##          : $create_index       => Create a BAM index when writing a coordinate-sorted BAM file
##          : $referencefile_path => Genome reference file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $referencefile_path;
    my $FILEHANDLE;

    ## Default(s)
    my $create_index;

    my $tmpl = {
        commands_ref => {
            default     => [],
            store       => \$commands_ref,
            strict_type => 1,
        },
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        FILEHANDLE => { store => \$FILEHANDLE, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = @{$commands_ref};

    if ( $create_index ne q{false} ) {

        push @commands, q{CREATE_INDEX=} . $create_index;
    }
    if ($referencefile_path) {

        push @commands, q{R=} . $referencefile_path;
    }
    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
