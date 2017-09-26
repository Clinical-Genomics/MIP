package MIP::Program::Variantcalling::Manta;

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
use File::Spec::Functions qw{ catdir catfile};

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
    our @EXPORT_OK = qw{ manta_config manta_workflow };

}

## Constants
Readonly my $SPACE => q{ };

sub manta_config {

    ## manta_config

    ## Function : Perl wrapper for writing Manta config recipe to $FILEHANDLE or return commands array. Based on Manta 1.0.0.
    ## Returns  : "@commands"
    ## Arguments: $infile_paths_ref, $referencefile_path, $outdirectory_path, $stderrfile_path, $stderrfile_path_append, $FILEHANDLE, $exome_analysis
    ##          : $infile_paths_ref       => Infile paths {REF}
    ##          : $referencefile_path     => Reference sequence file
    ##          : $outdirectory_path      => Outfile path
    ##          : $stderrfile_path        => Stderrfile path
    ##          : $stderrfile_path_append => Append stderr info to file path
    ##          : $FILEHANDLE             => Filehandle to write to
    ##          : $exome_analysis         => Set options for WES input: turn off depth filters

    my ($arg_href) = @_;

    ## Default(s)
    my $exome_analysis;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $referencefile_path;
    my $outdirectory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        outdirectory_path => { strict_type => 1, store => \$outdirectory_path },
        stderrfile_path   => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE     => { store => \$FILEHANDLE },
        exome_analysis => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$exome_analysis
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ configManta.py };

    #Reference sequence file
    if ($referencefile_path) {

        push @commands, q{--referenceFasta} . $SPACE . $referencefile_path;
    }

    if ($exome_analysis) {

        push @commands, q{--exome};
    }

    ## Infile
    push @commands, q{--bam} . $SPACE . join $SPACE . q{--bam} . $SPACE,
      @{$infile_paths_ref};

    if ($outdirectory_path) {

        push @commands, q{--runDir} . $SPACE . $outdirectory_path;
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

sub manta_workflow {

    ## manta_workflow

    ## Function : Perl wrapper for writing Manta workflow recipe to $FILEHANDLE or return commands array. Based on Manta 1.0.0.
    ## Returns  : "@commands"
    ## Arguments: $outdirectory_path, $stderrfile_path, $stderrfile_path_append, $outdirectory_path, $FILEHANDLE, $mode
    ##          : $outdirectory_path      => Outfile path
    ##          : $stderrfile_path        => Stderrfile path
    ##          : $stderrfile_path_append => Append stderr info to file path
    ##          : $FILEHANDLE             => Filehandle to write to
    ##          : $mode                   => Mode of parallel

    my ($arg_href) = @_;

    ## Default(s)
    my $mode;

    ## Flatten argument(s)
    my $outdirectory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        outdirectory_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outdirectory_path
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store => \$FILEHANDLE },
        mode       => {
            default     => q{local},
            allow       => [qw{ undef local sge }],
            strict_type => 1,
            store       => \$mode
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Manta
    my @commands = catfile( $outdirectory_path, q{runWorkflow.py} );

    ## Options
    if ($mode) {

        push @commands, q{--mode} . $SPACE . $mode;
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
