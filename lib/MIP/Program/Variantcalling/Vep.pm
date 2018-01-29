package MIP::Program::Variantcalling::Vep;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use strict;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ variant_effect_predictor variant_effect_predictor_install };
}

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

sub variant_effect_predictor {

## Function : Perl wrapper for writing variant_effect_predictor recipe to $FILEHANDLE or return commands array. Based on VEP 90.
## Returns  : @commands
## Arguments: $assembly                => Assembly version to use
##          : $buffer_size             => Sets the internal buffer size, corresponding to the number of variants that are read in to memory simultaneously
##          : $cache_directory         => VEP chache directory
##          : $distance                => Modify the distance up and/or downstream between a variant and a transcript for which VEP will assign the upstream_gene_variant or downstream_gene_variant consequences
##          : $FILEHANDLE              => Filehandle to write to
##          : $infile_path             => Infile path to read from
##          : $infile_format           => Input file format - one of "ensembl", "vcf", "hgvs", "id"
##          : $outfile_format          => Output file format
##          : $outfile_path            => Outfile path to write to
##          : $plugins_ref             => Use named plugin {REF}
##          : $reference_path          => Reference sequence file
##          : $regions_ref             => The regions to process {REF}
##          : $stderrfile_path         => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append  => Append stderr info to file path
##          : $stdoutfile_path         => Stdoutfile path
##          : $vep_features_ref        => Features to add to VEP

    my ($arg_href) = @_;

    ## Default(s)
    my $distance;
    my $fork;
    my $infile_format;
    my $outfile_format;

    ## Flatten argument(s)
    my $assembly;
    my $buffer_size;
    my $cache_directory;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $plugins_ref;
    my $reference_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $vep_features_ref;

    my $tmpl = {
        assembly    => { store => \$assembly, strict_type => 1, },
        buffer_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$buffer_size,
            strict_type => 1,
        },
        cache_directory => { store => \$cache_directory, strict_type => 1, },
        distance        => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 5000,
            store       => \$distance,
            strict_type => 1,
        },
        FILEHANDLE => { store => \$FILEHANDLE, },
        fork       => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$fork,
            strict_type => 1,
        },
        infile_format => {
            allow       => [qw{ ensembl vcf hgvs id }],
            default     => q{vcf},
            store       => \$infile_format,
            strict_type => 1,
        },
        infile_path    => { store => \$infile_path, strict_type => 1, },
        outfile_format => {
            allow       => [qw{ vcf tab json }],
            default     => q{vcf},
            store       => \$outfile_format,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        plugins_ref =>
          { default => [], store => \$plugins_ref, strict_type => 1, },
        reference_path => { store => \$reference_path, strict_type => 1, },
        regions_ref =>
          { default => [], store => \$regions_ref, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        vep_features_ref =>
          { default => [], store => \$vep_features_ref, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Vep
    # Stores commands depending on input parameters
    my @commands = q{vep};

    ## Options
    if ($fork) {

        push @commands, q{--fork} . $SPACE . $fork;
    }

    push @commands, q{--distance} . $SPACE . $distance;

    if ($buffer_size) {

        push @commands, q{--buffer_size} . $SPACE . $buffer_size;
    }
    if ($assembly) {

        push @commands, q{--assembly} . $SPACE . $assembly;
    }
    if ($reference_path) {

        push @commands, q{--fasta} . $SPACE . $reference_path;
    }
    if ($cache_directory) {

        push @commands, q{--dir_cache} . $SPACE . $cache_directory;
    }
    if ($infile_format) {

        push @commands, q{--format} . $SPACE . $infile_format;
    }
    if ($outfile_format) {

        push @commands, q{--} . $outfile_format;
    }

    # If regions limit output
    if ( @{$regions_ref} ) {

        push @commands, q{--chr} . $SPACE . join $COMMA, @{$regions_ref};
    }
    if ( @{$plugins_ref} ) {

        push @commands, q{--plugin} . $SPACE . join q{ --plugin },
          @{$plugins_ref};
    }
    if ( @{$vep_features_ref} ) {

        push @commands, q{--} . join q{ --}, @{$vep_features_ref};
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{--input_file} . $SPACE . $infile_path;
    }
    if ($outfile_path) {

        push @commands, q{--output_file} . $SPACE . $outfile_path;
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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub variant_effect_predictor_install {

## Function : Perl wrapper for vep INSTALL script. Based on version 90.
## Returns  : @commands

## Arguments: $assembly               => Assembly name to use if more than one during --AUTO
##          : $auto                   => Run installer without user prompts. Use "a" (API + Faidx/htslib),"l" (Faidx/htslib only), "c" (cache), "f" (FASTA), "p" (plugins) to specify parts to install.
##          : $cache_directory        => Set destination directory for cache files
##          : $FILEHANDLE             => Filehandle to write to
##          : $plugins_ref            => Vep plugins {REF}
##          : $species_ref            => Comma-separated list of species to install when using --AUTO {REF}
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $assembly;
    my $auto;
    my $cache_directory;
    my $FILEHANDLE;
    my $plugins_ref;
    my $species_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        assembly        => { store => \$assembly,        strict_type => 1, },
        auto            => { store => \$auto,            strict_type => 1, },
        cache_directory => { store => \$cache_directory, strict_type => 1, },
        FILEHANDLE      => { store => \$FILEHANDLE, },
        plugins_ref =>
          { default => [], store => \$plugins_ref, strict_type => 1, },
        species_ref => {
            default     => [qw{ homo_sapiens }],
            store       => \$species_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ perl INSTALL.pl };

    if ($auto) {

        push @commands, q{--AUTO} . $SPACE . $auto;
    }
    if ( @{$plugins_ref} ) {

        push @commands, q{--PLUGINS} . $SPACE . join $COMMA, @{$plugins_ref};
    }
    if ($cache_directory) {

        push @commands, q{--CACHEDIR} . $SPACE . $cache_directory;
    }
    if ( @{$species_ref} ) {

        push @commands, q{--SPECIES} . $SPACE . join $COMMA, @{$species_ref};
    }
    if ($assembly) {

        push @commands, q{--ASSEMBLY} . $SPACE . $assembly;
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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,
        }
    );

    return @commands;
}

1;
