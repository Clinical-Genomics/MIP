package MIP::Program::Multiqc;

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

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ multiqc };
}

sub multiqc {

## Function  : Perl wrapper for writing multiqc recipe to already open $filehandle or return commands array. Based on multiqc 0.8.dev0.
## Returns   : @commands
## Arguments : $filehandle             => Filehandle to write to
##           : $force                  => Force overwrite of output files
##           : $indir_path             => Indir path
##           : $outdir_path            => Outdir path
##           : $stderrfile_path        => Stderrfile path
##           : $stderrfile_path_append => Append stderr info to file
##           : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $indir_path;
    my $outdir_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $force;

    my $tmpl = {
        filehandle => { store => \$filehandle, },
        force      => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$force,
            strict_type => 1,
        },
        indir_path => {
            defined     => 1,
            required    => 1,
            store       => \$indir_path,
            strict_type => 1,
        },
        outdir_path     => { store => \$outdir_path,     strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => q{multiqc}, } ), );

    if ($force) {

        push @commands, q{--force};
    }

    if ($outdir_path) {

        push @commands, q{--outdir} . $SPACE . $outdir_path;
    }

    push @commands, $indir_path;

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
