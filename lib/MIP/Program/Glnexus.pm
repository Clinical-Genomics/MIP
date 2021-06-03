package MIP::Program::Glnexus;

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

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ glnexus_merge };
}

sub glnexus_merge {

## Function : Perl wrapper for generic commands module
## Returns  : @commands
## Arguments: $config                 => Allows us to process vcf files from different variant callers.
##          : $dir                    => Directory to write temporary files into
##          : $infile_paths_ref       => Space separated list of input vcf files
##          : $memory                 => Memory limit
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $threads => Number of threads to use
    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $config;
    my $dir;
    my $filehandle;
    my $infile_paths_ref;
    my $memory;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $threads;

    my $tmpl = {
        config => {
            allow => [
                qw{ gatk gatk_unfiltered xAtlas xAtlas_unfiltered weCall weCall_unfiltered
                  DeepVariant DeepVariantWGS DeepVariantWES DeepVariant_unfiltered Strelka2 }
            ],
            required    => 1,
            store       => \$config,
            strict_type => 1,
        },
        dir => {
            required    => 1,
            store       => \$dir,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        memory => {
            store       => \$memory,
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
        threads => {
            store       => \$threads,
            strict_type => 1,
        }
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => q{glnexus_cli}, } ), );

    if ($threads) {
        push @commands, q{--threads} . $SPACE . $threads;
    }
    if ($memory) {
        push @commands, q{--mem-gbytes} . $SPACE . $memory;
    }

    push @commands, q{--config} . $SPACE . $config;
    push @commands, q{--dir} . $SPACE . $dir;
    push @commands, join $SPACE, @{$infile_paths_ref};

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
