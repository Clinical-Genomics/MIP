package MIP::Program::Variant_integrity;

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
use MIP::Constants qw{ $DASH $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ variant_integrity_father variant_integrity_mendel };

}

sub variant_integrity_father {

## Function : Perl wrapper for writing Variant_integrity father recipe to $filehandle or return commands array. Based on Variant_integrity 3.7.0.
## Returns  : @commands
## Arguments: $case_file              => Family file
##          : $case_type              => Setup of family file
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $verbosity              => Increase output verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_file;
    my $case_type;
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $verbosity;

    my $tmpl = {
        case_file => {
            defined     => 1,
            required    => 1,
            store       => \$case_file,
            strict_type => 1,
        },
        case_type => {
            allow       => [qw{ ped alt cmms mip }],
            store       => \$case_type,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path,    strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        verbosity => {
            allow       => qr{ \A \w+ \z }sxm,
            store       => \$verbosity,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ variant_integrity };

    if ($verbosity) {

        push @commands, $DASH . $verbosity;
    }

    if ($case_file) {

        push @commands, q{--family_file} . $SPACE . $case_file;
    }

    if ($case_type) {

        push @commands, q{--family_type} . $SPACE . $case_type;
    }

    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub variant_integrity_mendel {

## Function : Perl wrapper for writing Variant_integrity mendel recipe to $filehandle or return commands array. Based on Variant_integrity 3.7.0.
## Returns  : @commands
## Arguments: $case_file              => Family file
##          : $case_type              => Setup of family file
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $verbosity              => Increase output verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_file;
    my $case_type;
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $verbosity;

    my $tmpl = {
        case_file => {
            defined     => 1,
            required    => 1,
            store       => \$case_file,
            strict_type => 1,
        },
        case_type => {
            allow       => [qw{ ped alt cmms mip }],
            store       => \$case_type,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path,    strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        verbosity => {
            allow       => qr{ \A \w+ \z }sxm,
            store       => \$verbosity,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ variant_integrity };

    if ($verbosity) {

        push @commands, $DASH . $verbosity;
    }

    if ($case_file) {

        push @commands, q{--family_file} . $SPACE . $case_file;
    }

    if ($case_type) {

        push @commands, q{--family_type} . $SPACE . $case_type;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
