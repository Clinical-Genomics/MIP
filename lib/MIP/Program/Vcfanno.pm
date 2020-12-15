package MIP::Program::Vcfanno;

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
    our @EXPORT_OK = qw{ vcfanno };
}

sub vcfanno {

## Function : Perl wrapper for writing vcfanno recipe to $filehandle or return commands array. Based on vcfanno 0.1.0.
## Returns  : @commands
## Arguments: $ends                   => annotate the start and end as well as the interval itself
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $luafile_path           => Optional path to a file containing custom javascript functions to be used as ops
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $toml_configfile_path   => Toml config file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $ends;
    my $filehandle;
    my $infile_path;
    my $luafile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $toml_configfile_path;

    my $tmpl = {
        ends => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$ends,
            strict_type => 1,
        },
        filehandle   => { store => \$filehandle, },
        infile_path  => { store => \$infile_path, strict_type => 1, },
        luafile_path => { store => \$luafile_path, strict_type => 1, },
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
        toml_configfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$toml_configfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => q{vcfanno}, } ), );

    if ($luafile_path) {

        push @commands, q{-lua} . $SPACE . $luafile_path;
    }
    if ($ends) {

        push @commands, q{-ends};
    }

    push @commands, $toml_configfile_path;

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

1;
