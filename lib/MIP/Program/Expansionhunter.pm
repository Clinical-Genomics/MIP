package MIP::Program::Expansionhunter;

use 5.026;
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
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ expansionhunter };
}

sub expansionhunter {

## Function : Perl wrapper for writing Expansion Hunter command to $filehandle or return commands array. Based on Expansionhunter version 3.0.0
## Returns  : @commands
## Arguments: $filehandle                => Filehandle to write to
##          : $infile_path               => Path to sorted, indexed bam file
##          : $log_level                 => Log level
##          : $outfile_path_prefix       => Path to vcf file
##          : $reference_genome_path     => Path to reference fasta file
##          : $region_extension_length   => Specifies how far from on/off-target regions to search for informative reads
##          : $sex                       => Sex of the sample (Defaults to female)
##          : $stderrfile_path           => Stderrfile path
##          : $stderrfile_path_append    => Append stderr info to file path
##          : $stdoutfile_path           => Stdoutfile path
##          : $variant_catalog_file_path => Path to dir with repeat-specification files

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $log_level;
    my $outfile_path_prefix;
    my $reference_genome_path;
    my $region_extension_length;
    my $sex;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $variant_catalog_file_path;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        log_level => {
            allow       => [qw{ trace debug info warn error }],
            default     => q{info},
            store       => \$log_level,
            strict_type => 1,
        },
        outfile_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path_prefix,
            strict_type => 1,
        },
        reference_genome_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_path,
            strict_type => 1,
        },
        region_extension_length => {
            allow       => qr/ ^\d+$ /xms,
            defined     => 1,
            store       => \$region_extension_length,
            strict_type => 1,
        },
        sex => {
            allow       => [ undef, qw{ female male } ],
            store       => \$sex,
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
        variant_catalog_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$variant_catalog_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{ExpansionHunter};

    ## Required arguments
    push @commands, q{--reads} . $SPACE . $infile_path;

    push @commands, q{--log-level} . $SPACE . $log_level;

    push @commands, q{--reference} . $SPACE . $reference_genome_path;

    push @commands, q{--variant-catalog} . $SPACE . $variant_catalog_file_path;

    push @commands, q{--output-prefix} . $SPACE . $outfile_path_prefix;

    if ($region_extension_length) {

        push @commands, q{--region-extension-length} . $SPACE . $region_extension_length;
    }

    if ($sex) {

        push @commands, q{--sex} . $SPACE . $sex;
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
