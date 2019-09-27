package MIP::Program::Variantcalling::Vep;

use Carp;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use Params::Check qw{ check allow last_error };
use charnames qw{ :full :short };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPANM
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ variant_effect_predictor variant_effect_predictor_install };
}

## Constants
Readonly my $LENGTH_CHR_1 => 248_956_422;

sub variant_effect_predictor {

## Function : Perl wrapper for writing variant_effect_predictor recipe to $FILEHANDLE or return commands array. Based on VEP 90.
## Returns  : @commands
## Arguments: $assembly                => Assembly version to use
##          : $buffer_size             => Sets the internal buffer size, corresponding to the number of variants that are read in to memory simultaneously
##          : $cache_directory         => VEP chache directory
##          : $custom_annotations_ref  => Custom annotations {REF}
##          : $distance                => Modify the distance up and/or downstream between a variant and a transcript for which VEP will assign the upstream_gene_variant or downstream_gene_variant consequences
##          : $FILEHANDLE              => Filehandle to write to
##          : $infile_path             => Infile path to read from
##          : $infile_format           => Input file format - one of "ensembl", "vcf", "hgvs", "id"
##          : $max_sv_size             => Extend the maximum Structural Variant size VEP can process
##          : $outfile_format          => Output file format
##          : $outfile_path            => Outfile path to write to
##          : $plugins_dir_path        => Path to plugins directory
##          : $plugins_ref             => Use named plugin {REF}
##          : $reference_path          => Reference sequence file
##          : $regions_ref             => The regions to process {REF}
##          : $stderrfile_path         => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append  => Append stderr info to file path
##          : $stdoutfile_path         => Stdoutfile path
##          : $synonyms_file_path      => Contig synonyms
##          : $vep_features_ref        => Features to add to VEP

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $assembly;
    my $buffer_size;
    my $cache_directory;
    my $custom_annotations_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $plugins_dir_path;
    my $plugins_ref;
    my $reference_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $synonyms_file_path;
    my $vep_features_ref;

    ## Default(s)
    my $distance;
    my $fork;
    my $infile_format;
    my $max_sv_size;
    my $outfile_format;

    my $tmpl = {
        assembly => {
            store       => \$assembly,
            strict_type => 1,
        },
        buffer_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$buffer_size,
            strict_type => 1,
        },
        cache_directory => {
            store       => \$cache_directory,
            strict_type => 1,
        },
        custom_annotations_ref => {
            default     => [],
            store       => \$custom_annotations_ref,
            strict_type => 1,
        },
        distance => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 5000,
            store       => \$distance,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        fork => {
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
        infile_path => {
            store       => \$infile_path,
            strict_type => 1,
        },
        max_sv_size => {
            allow       => qr/\A \d+ \z/xsm,
            default     => $LENGTH_CHR_1,
            store       => \$max_sv_size,
            strict_type => 1,
        },
        outfile_format => {
            allow       => [qw{ vcf tab json }],
            default     => q{vcf},
            store       => \$outfile_format,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        plugins_dir_path => {
            store       => \$plugins_dir_path,
            strict_type => 1,
        },
        plugins_ref => {
            default     => [],
            store       => \$plugins_ref,
            strict_type => 1,
        },
        reference_path => {
            store       => \$reference_path,
            strict_type => 1,
        },
        regions_ref => {
            default     => [],
            store       => \$regions_ref,
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
        synonyms_file_path => {
            store       => \$synonyms_file_path,
            strict_type => 1,
        },
        vep_features_ref => {
            default     => [],
            store       => \$vep_features_ref,
            strict_type => 1,
        },
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

    push @commands, q{--max_sv_size} . $SPACE . $max_sv_size;

    if ($outfile_format) {

        push @commands, q{--} . $outfile_format;
    }

    # If regions limit output
    if ( @{$regions_ref} ) {

        push @commands, q{--chr} . $SPACE . join $COMMA, @{$regions_ref};
    }
    if ($synonyms_file_path) {

        push @commands, q{--synonyms} . $SPACE . $synonyms_file_path;
    }
    if ($plugins_dir_path) {
        push @commands, q{--dir_plugins} . $SPACE . $plugins_dir_path;
    }
    if ( @{$plugins_ref} ) {

        push @commands, q{--plugin} . $SPACE . join q{ --plugin }, @{$plugins_ref};
    }
    if ( @{$custom_annotations_ref} ) {

        push @commands, q{--custom} . $SPACE . join q{ --custom },
          @{$custom_annotations_ref};
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
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub variant_effect_predictor_install {

## Function : Perl wrapper for vep INSTALL script. Based on version 92.
## Returns  : @commands
## Arguments: $assembly               => Assembly name to use if more than one during --AUTO
##          : $auto                   => Run installer without user prompts. Use "a" (API + Faidx/htslib),"l" (Faidx/htslib only), "c" (cache), "f" (FASTA), "p" (plugins) to specify parts to install.
##          : $cache_directory        => Set destination directory for cache files
##          : $cache_version          => Set cache version to download
##          : $FILEHANDLE             => Filehandle to write to
##          : $no_update              => Don't update
##          : $no_htslib              => Don't attempt to install Bio::DB::HTS/htslib
##          : $plugins_ref            => Vep plugins {REF}
##          : $species_ref            => Comma-separated list of species to install when using --AUTO {REF}
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $version                => Version to install

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $assembly;
    my $auto;
    my $cache_directory;
    my $cache_version;
    my $FILEHANDLE;
    my $no_htslib;
    my $no_update;
    my $plugins_ref;
    my $species_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $version;

    my $tmpl = {
        assembly        => { store => \$assembly,        strict_type => 1, },
        auto            => { store => \$auto,            strict_type => 1, },
        cache_directory => { store => \$cache_directory, strict_type => 1, },
        cache_version   => { store => \$cache_version,   strict_type => 1, },
        FILEHANDLE      => { store => \$FILEHANDLE, },
        no_update       => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            store       => \$no_update,
            strict_type => 1,
        },
        no_htslib => {
            allow       => [ undef, 0, 1 ],
            store       => \$no_htslib,
            strict_type => 1,
        },
        plugins_ref => { default => [], store => \$plugins_ref, strict_type => 1, },
        species_ref => {
            default     => [qw{ homo_sapiens }],
            store       => \$species_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        version         => { store => \$version,         strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ INSTALL.pl };

    if ($auto) {

        push @commands, q{--AUTO} . $SPACE . $auto;
    }
    if ($version) {

        push @commands, q{--VERSION} . $SPACE . $version;
    }
    if ($cache_version) {

        push @commands, q{--CACHE_VERSION} . $SPACE . $cache_version;
    }
    if ($no_update) {

        push @commands, q{--NO_UPDATE};
    }
    if ($no_htslib) {

        push @commands, q{--NO_HTSLIB};
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
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

1;
