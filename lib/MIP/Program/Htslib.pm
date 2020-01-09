package MIP::Program::Htslib;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ htslib_bgzip htslib_tabix };
}

sub htslib_bgzip {

## Function : Perl wrapper for writing bgzip recipe to $filehandle or return commands array. Based on htslib 1.3.1.
## Returns  : @commands
## Arguments: $decompress             => Decompress file
##          : $force                  => Force
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $write_to_stdout        => Write on standard output, keep original files unchanged

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $decompress;
    my $force;
    my $write_to_stdout;

    my $tmpl = {
        decompress => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$decompress
        },
        force => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$force
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path     => { strict_type => 1, store => \$infile_path },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
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

    if ($decompress) {

        push @commands, q{--decompress};
    }
    if ($force) {

        push @commands, q{--force};
    }

    if ($write_to_stdout) {

        push @commands, q{--stdout};
    }

    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
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

sub htslib_tabix {

## Function : Perl wrapper for writing tabix recipe to $filehandle or return commands array. Based on htslib 1.3.1.
## Returns  : @commands
## Arguments: $begin                  => Column number for region start
##          : $end                    => Column number for region end
##          : $filehandle             => Filehandle to write to
##          : $force                  => Overwrite existing index without asking
##          : $infile_path            => Infile path to read from
##          : $preset                 => Preset
##          : $regions_ref            => Regions to process {REF}
##          : $sequence               => Column number for sequence names
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $with_header            => Include header
##          : $zero_based             => Coordinates are zero-based

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $begin;
    my $end;
    my $force;
    my $preset;
    my $sequence;
    my $with_header;
    my $zero_based;

    my $tmpl = {
        begin => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 0,
            store       => \$begin,
            strict_type => 1,
        },
        end => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 0,
            store       => \$end,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        force => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$force,
            strict_type => 1,
        },
        infile_path => { store => \$infile_path, strict_type => 1, },
        preset      => {
            allow       => [ undef, qw{ gff bed sam vcf } ],
            store       => \$preset,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        sequence    => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 0,
            store       => \$sequence,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        with_header => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$with_header,
            strict_type => 1,
        },
        zero_based => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$zero_based,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{tabix};

    if ($force) {

        push @commands, q{--force};
    }

    if ($preset) {

        push @commands, q{--preset} . $SPACE . $preset;
    }

    if ($with_header) {

        push @commands, q{--print-header};
    }

    if ($zero_based) {

        push @commands, q{--zero-based};
    }
    if ($begin) {

        push @commands, q{--begin} . $SPACE . $begin;
    }

    if ($end) {

        push @commands, q{--end} . $SPACE . $end;
    }

    if ($sequence) {

        push @commands, q{--sequence} . $SPACE . $sequence;
    }

    if ($infile_path) {

        push @commands, $infile_path;
    }

    if ( @{$regions_ref} ) {

        push @commands, join $SPACE, @{$regions_ref};
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
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
