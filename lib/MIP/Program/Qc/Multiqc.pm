package MIP::Program::Qc::Multiqc;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ multiqc };
}

## Constants
Readonly my $SPACE => q{ };

sub multiqc {

## Function  : Perl wrapper for writing multiqc recipe to already open $filehandle or return commands array. Based on multiqc 0.8.dev0.
## Returns   : "@commands"
## Arguments : $indir_path             => Indir path
##           : $outdir_path            => Outdir path
##           : $stderrfile_path        => Stderrfile path
##           : $stdoutfile_path        => Stdoutfile path
##           : $filehandle             => Filehandle to write to
##           : $stderrfile_path_append => Append stderr info to file
##           : $force                  => Force overwrite of output files

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $indir_path;
    my $outdir_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $filehandle;

    ## Default(s)
    my $force;

    my $tmpl = {
        indir_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$indir_path
        },
        outdir_path     => { strict_type => 1, store => \$outdir_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        filehandle => { store => \$filehandle },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        force => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$force
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{multiqc};

    ## Options
    if ($force) {

        push @commands, q{--force};
    }

    ## Outdir
    if ($outdir_path) {

        push @commands, q{--outdir} . $SPACE . $outdir_path;
    }

    ## Indir
    push @commands, $indir_path;

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

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
