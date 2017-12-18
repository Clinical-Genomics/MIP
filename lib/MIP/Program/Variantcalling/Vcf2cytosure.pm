package MIP::Program::Variantcalling::Vcf2cytosure;

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
use autodie qw{ :all };
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
    our @EXPORT_OK = qw{ vcf2cytosure_convert };
}

## Constants
Readonly my $SPACE => q{ };

sub vcf2cytosure_convert {

## Function : Perl wrapper for Vcf2cytosure 0.2.0.
## Returns  : @commands
## Arguments: $coverage_file          => Path to coverage file
##          : $FILEHANDLE             => Filehandle to write to
##          : $frequency              => Maximum frequency
##          : $frequency_tag          => Frequency tag of the info field
##          : $infile_paths_ref       => VCF infiles paths {REF}
##          : $no_filter              => Disable any filtering
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $variant_size           => Minimum variant size.
##          : $vcf_infile_path        => Path to VCF infile

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $coverage_file;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $vcf_infile_path;

    ## Default(s)
    my $frequency;
    my $frequency_tag;
    my $no_filter;
    my $variant_size;

    my $tmpl = {
        coverage_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$coverage_file,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        frequency => {
            default     => 0.01,
            strict_type => 1,
            store       => \$frequency,
        },
        frequency_tag => {
            default     => q{FRQ},
            strict_type => 1,
            store       => \$frequency_tag,
        },
        no_filter => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_filter,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
        variant_size => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 5000,
            store       => \$variant_size,
            strict_type => 1,
        },
        vcf_infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$vcf_infile_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{vcf2cytosure};

    # Option: specify minimum variant size
    if ($variant_size) {
        push @commands, q{--size} . $SPACE . $variant_size;
    }

    # Option: specify maximum frequency
    if ($frequency) {
        push @commands, q{--frequency} . $SPACE . $frequency;
    }

    # Option: specify frequency tag
    if ($frequency_tag) {
        push @commands, q{--frequency_tag} . $SPACE . $frequency_tag;
    }

    # Option: no filtering. Overrides previous filtering options
    if ($no_filter) {
        push @commands, q{--no-filter};
    }

    # Coverage file
    push @commands, q{--coverage} . $SPACE . $coverage_file;

    # VCF file
    push @commands, q{--vcf} . $SPACE . $vcf_infile_path;

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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
