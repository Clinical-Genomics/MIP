package MIP::Program::Variantcalling::Variant_integrity;

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

use FindBin qw{ $Bin };    #Find directory of script
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ variant_integrity_mendel variant_integrity_father };

}

## Constants
Readonly my $SPACE => q{ };
Readonly my $DASH  => q{-};

sub variant_integrity_mendel {

## variant_integrity_mendel

## Function : Perl wrapper for writing Variant_integrity mendel recipe to $FILEHANDLE or return commands array. Based on Variant_integrity 3.7.0.
## Returns  : "@commands"

## Arguments: $infile_path, $family_file, $outfile_path, $stderrfile_path, $stderrfile_path_append, $FILEHANDLE, $verbosity, $family_type
##          : $infile_path            => Infile path to read from
##          : $family_file            => Family file
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $verbosity              => Increase output verbosity
##          : $family_type            => Setup of family file

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;

    ## Flatten argument(s)
    my $infile_path;
    my $family_file;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $family_type;

    my $tmpl = {

        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        family_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_file
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE  => { store => \$FILEHANDLE },
        family_type => {
            allow       => [qw{ ped alt cmms mip }],
            strict_type => 1,
            store       => \$family_type
        },
        verbosity => {
            allow       => qr/^\w+$/,
            strict_type => 1,
            store       => \$verbosity
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ variant_integrity };

    ## Options
    if ($verbosity) {

        push @commands, $DASH . $verbosity;
    }

    if ($family_file) {

        push @commands, q{--family_file} . $SPACE . $family_file;
    }

    if ($family_type) {

        push @commands, q{--family_type} . $SPACE . $family_type;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands, q{mendel};

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub variant_integrity_father {

## variant_integrity_father

## Function : Perl wrapper for writing Variant_integrity father recipe to $FILEHANDLE or return commands array. Based on Variant_integrity 3.7.0.
## Returns  : "@commands"
## Arguments: $infile_path, $family_file, $outfile_path, $stderrfile_path, $stderrfile_path_append, $FILEHANDLE, $verbosity, $family_type
##          : $infile_path            => Infile path to read from
##          : $family_file            => Family file
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $verbosity              => Increase output verbosity
##          : $family_type            => Setup of family file

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;

    ## Flatten argument(s)
    my $infile_path;
    my $family_file;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $family_type;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        family_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$family_file
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE  => { store => \$FILEHANDLE },
        family_type => {
            allow       => [qw{ ped alt cmms mip }],
            strict_type => 1,
            store       => \$family_type
        },
        verbosity => {
            allow       => qr/^\w+$/,
            strict_type => 1,
            store       => \$verbosity
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ variant_integrity };

    ## Options
    if ($verbosity) {

        push @commands, $DASH . $verbosity;
    }

    if ($family_file) {

        push @commands, q{--family_file} . $SPACE . $family_file;
    }

    if ($family_type) {

        push @commands, q{--family_type} . $SPACE . $family_type;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands, q{father};

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

1;
