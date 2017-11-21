package MIP::Program::Utility::Htslib;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ htslib_bgzip htslib_tabix };
}

## Constants
Readonly my $SPACE => q{ };

sub htslib_bgzip {

## Function : Perl wrapper for writing bgzip recipe to $FILEHANDLE or return commands array. Based on htslib 1.3.1.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $infile_path            => Infile path to read from
##          : $decompress             => Decompress file
##          : $write_to_stdout        => Write on standard output, keep original files unchanged

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $infile_path;
    my $outfile_path;

    ## Default(s)
    my $decompress;
    my $write_to_stdout;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        infile_path  => { strict_type => 1, store => \$infile_path },
        outfile_path => { strict_type => 1, store => \$outfile_path },
        decompress   => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$decompress
        },
        write_to_stdout => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$write_to_stdout
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{bgzip};

    ## Options
    if ($decompress) {

        push @commands, q{--decompress};
    }

    if ($write_to_stdout) {

        push @commands, q{--stdout};
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    #Specify output filename
    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }

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

sub htslib_tabix {

## Function : Perl wrapper for writing tabix recipe to $FILEHANDLE or return commands array. Based on htslib 1.3.1.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $regions_ref            => The regions to process {REF}
##          : $force                  => Overwrite existing index without asking
##          : $preset                 => Preset
##          : $with_header            => Include header

    my ($arg_href) = @_;

    ## Default(s)
    my $with_header;
    my $force;
    my $preset;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path  => { strict_type => 1, store => \$infile_path },
        outfile_path => { strict_type => 1, store => \$outfile_path },
        force        => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$force
        },
        preset => {
            allow       => [ undef, qw{ gff bed sam vcf } ],
            strict_type => 1,
            store       => \$preset
        },
        with_header => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$with_header
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{tabix};

    ## Options
    if ($force) {

        push @commands, q{--force};
    }

    if ($preset) {

        push @commands, q{--preset} . $SPACE . $preset;
    }

    #Include header
    if ($with_header) {

        push @commands, q{--print-header};
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    #Limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, join $SPACE, @{$regions_ref};
    }

    #Specify output filename
    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }

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
