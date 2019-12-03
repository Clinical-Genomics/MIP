package MIP::Program::Manta;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile};
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ manta_config manta_workflow };

}

sub manta_config {

## Function : Perl wrapper for writing Manta config recipe to $filehandle or return commands array. Based on Manta 1.5.0.
## Returns  : @commands
## Arguments: $call_regions_file_path => Call regions file path
##          : $exome_analysis         => Set options for WES input: turn off depth filters
##          : $filehandle             => Filehandle to write to
##          : $infile_paths_ref       => Infile paths {REF}
##          : $outdirectory_path      => Outfile path
##          : $referencefile_path     => Reference sequence file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $call_regions_file_path;
    my $filehandle;
    my $infile_paths_ref;
    my $outdirectory_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $exome_analysis;

    my $tmpl = {
        call_regions_file_path => {
            store       => \$call_regions_file_path,
            strict_type => 1,
        },
        exome_analysis => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$exome_analysis,
            strict_type => 1,
        },
        filehandle       => { store => \$filehandle, },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outdirectory_path  => { store => \$outdirectory_path, strict_type => 1, },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ configManta.py };

    if ($referencefile_path) {

        push @commands, q{--referenceFasta} . $SPACE . $referencefile_path;
    }

    if ($call_regions_file_path) {

        push @commands, q{--callRegions} . $SPACE . $call_regions_file_path;
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub manta_workflow {

## Function : Perl wrapper for writing Manta workflow recipe to $filehandle or return commands array. Based on Manta 1.5.0.
## Returns  : "@commands"
## Arguments: $core_number            => Number of cores to use
##          : $filehandle             => Filehandle to write to
##          : $mode                   => Mode of parallel
##          : $outdirectory_path      => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $core_number;
    my $filehandle;
    my $outdirectory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $mode;

    my $tmpl = {
        core_number => {
            allow       => qr{ \A\d+\z }sxm,
            store       => \$core_number,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        mode       => {
            allow       => [qw{ undef local sge }],
            default     => q{local},
            store       => \$mode,
            strict_type => 1,
        },
        outdirectory_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdirectory_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Manta
    my @commands = qw{ runWorkflow.py };

    push @commands, catfile( $outdirectory_path, q{runWorkflow.py} );

    if ($mode) {

        push @commands, q{--mode} . $SPACE . $mode;
    }

    if ($core_number) {

        push @commands, q{-j} . $SPACE . $core_number;
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
