package MIP::Program::Qc::Fastqc;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(fastqc);

}

## Constants
Readonly my $SPACE => q{ };

sub fastqc {

##fastqc

##Function : Perl wrapper for writing fastqc recipe to already open $filehandle or return commands array. Based on cp 0.11.5
##Returns  : "@commands"
##Arguments: $filehandle, $infile_path, $outdirectory_path, $extract, $quiet
##         : $filehandle        => Filehandle to write to
##         : $infile_path       => Infile path
##         : $outdirectory_path => Outdirectory path
##         : $extract           => If set then the zipped output file will be uncompressed in the same directory after it has been created
##         : $quiet             => Supress all progress messages on stdout and only report errors

    my ($arg_href) = @_;

    ## Default(s)
    my $extract;
    my $quiet;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outdirectory_path;

    my $tmpl = {
        filehandle  => { store => \$filehandle },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outdirectory_path => { strict_type => 1, store => \$outdirectory_path },
        extract           => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$extract
        },
        quiet => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Fastqc
    # Stores commands depending on input parameters
    my @commands = qw(fastqc);

    if ($quiet) {

        # Supress all progress messages on stdout and only report errors
        push @commands, q{--quiet};
    }
    if ($extract) {

# Zipped output file will be uncompressed in the same directory after it has been created
        push @commands, q{--extract};
    }
    if ($outdirectory_path) {

        # Create all output files in the specified output directory
        push @commands, q{--outdir} . $SPACE . $outdirectory_path;
    }

    ## Infile
    push @commands, $infile_path;

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            filehandle   => $filehandle,
        }
    );
    return @commands;
}

1;
