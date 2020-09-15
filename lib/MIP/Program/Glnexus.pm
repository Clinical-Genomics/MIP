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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ glnexus_merge };
}

sub glnexus_merge {

## Function : Perl wrapper for generic commands module
## Returns  : @commands
## Arguments: $config                 => Allows us to process vcf files from different variant callers.
##          : $infile_paths_ref            => Space separated list of input vcf files
##          : $stdoutfile_path        => Stdoutfile path
    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $config;
    my $filehandle;
    my $infile_paths_ref;
    my $stdoutfile_path;

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
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ glnexus_cli };

    push @commands, q{--config} . $SPACE . $config;
    push @commands, join $SPACE, @{$infile_paths_ref};

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
