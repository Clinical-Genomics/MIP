package MIP::Program::Download::Download_reference;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile canonpath };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_reference };
}

## Constants
Readonly my $SPACE => q{ };

sub download_reference {

## Function : Perl wrapper for writing download_reference command to filehandle.
##          : Based on download_reference version 0.0.4
## Returns  : @commands
## Arguments: $config_file_path              => Config file path
##          : $filehandle                    => Filehandle to write to
##          : $pipeline                      => Pipeline
##          : $reference_dir_path            => Reference directory
##          : $reference_genome_versions_ref => Array with genome versions to downlaod {REF}
##          : $stderrfile_path               => Stderrfile path
##          : $stderrfile_path_append        => Append stderr info to file path
##          : $stdoutfile_path               => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $config_file_path;
    my $filehandle;
    my $pipeline;
    my $reference_dir_path;
    my $reference_genome_versions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        config_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$config_file_path,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        pipeline => {
            allow       => [qw{ rare_disease rna }],
            defined     => 1,
            required    => 1,
            store       => \$pipeline,
            strict_type => 1,
        },
        reference_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir_path,
            strict_type => 1,
        },
        reference_genome_versions_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_versions_ref,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{mip download};

    push @commands, $pipeline;

    push @commands, q{--config_file} . $SPACE . $config_file_path;

    push @commands,
      q{--reference_dir} . $SPACE . canonpath($reference_dir_path);

    push @commands,
        q{--reference_genome_versions}
      . $SPACE
      . join $SPACE
      . q{--reference_genome_versions}
      . $SPACE,
      @{$reference_genome_versions_ref};

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
