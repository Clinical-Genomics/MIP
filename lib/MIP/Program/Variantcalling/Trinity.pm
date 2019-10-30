package MIP::Program::Variantcalling::Trinity;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.0;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ trinity_genome_guided };

}

## Constants
Readonly my $SPACE => q{ };

sub trinity_genome_guided {

## Function : Perl wrapper for writing Trinity genome guided asssembly recipe to $filehandle or return commands array.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $max_intron_distance    => Maximum intron distance
##          : $max_memory             => Maximum memory (gigabytes)
##          : $number_cores           => Number of CPU cores
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $max_intron_distance;
    my $max_memory;
    my $number_cpu;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle  => { store => \$filehandle },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        max_intron_distance => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 10_000,
            store       => \$max_intron_distance,
            strict_type => 1,
        },
        max_memory => {
            default     => 50,
            store       => \$max_memory,
            strict_type => 1,
        },
        number_cpu =>
          { default => 16, store => \$number_cpu, strict_type => 1, },
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Trinity
    my @commands = qw{ Trinity };

    ## Infile
    push @commands, q{--genome_guided_bam} . $SPACE . $infile_path;

    ## cpu
    push @commands,
      q{--genome_guided_max_intron} . $SPACE . $max_intron_distance;

    ## intron distance
    push @commands, q{--CPU} . $SPACE . $number_cpu;

    ## maximum memory
    push @commands, q{--max_memory} . $SPACE . $max_memory . q{G};

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    push @commands,
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
