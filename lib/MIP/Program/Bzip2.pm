package MIP::Program::Bzip2;

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
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ bzip2 };
}

sub bzip2 {

## Function : Perl wrapper for writing bzip2 recipe to $filehandle or return commands array. Based on bzip2 v1.0.6
## Returns  : @commands
## Arguments: $decompress             => Decompress bzip2 file
##          : $filehandle             => Filehandle to write to
##          : $force                  => Overwrite of output files
##          : $infile_path            => Infile path
##          : $outfile_path           => Path to output file
##          : $quiet                  => Suppress all warnings
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdout                 => Write on standard output
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $decompress;
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $force;
    my $quiet;
    my $stdout;
    my $verbose;

    my $tmpl = {
        decompress => {
            store       => \$decompress,
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
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$quiet,
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
        stdout => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$stdout,
            strict_type => 1,
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => q{bzip2}, } ), );

    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($verbose) {
        push @commands, q{--verbose};
    }

    if ($force) {
        push @commands, q{--force};
    }

    if ($decompress) {
        push @commands, q{--decompress};
    }

    ## Write to stdout stream
    if ($stdout) {
        push @commands, q{--stdout};
    }

    ## Infile
    push @commands, $infile_path;

    ## Outfile
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
