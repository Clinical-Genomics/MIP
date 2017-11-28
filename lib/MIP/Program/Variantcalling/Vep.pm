package MIP::Program::Variantcalling::Vep;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };

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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ variant_effect_predictor variant_effect_predictor_install };
}

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

sub variant_effect_predictor {

##variant_effect_predictor

##Function : Perl wrapper for writing variant_effect_predictor recipe to $FILEHANDLE or return commands array. Based on VEP 90.
##Returns  : @commands
##Arguments: $plugins_ref, $regions_ref, $outfile_path, $infile_path, $stdoutfile_path, $stderrfile_path, $stderrfile_path_append, $FILEHANDLE, $reference_path, $script_path, $assembly, $buffer_size, $infile_format, $outfile_format, $fork, $no_progress, $offline
##         : $plugins_ref             => Use named plugin {REF}
##         : $regions_ref             => The regions to process {REF}
##         : $vep_features_ref        => Features to add to VEP
##         : $outfile_path            => Outfile path to write to
##         : $infile_path             => Infile path to read from
##         : $stdoutfile_path         => Stdoutfile path
##         : $stderrfile_path         => Stderr file path to write to {OPTIONAL}
##         : $stderrfile_path_append  => Append stderr info to file path
##         : $FILEHANDLE              => Filehandle to write to
##         : $reference_path          => Reference sequence file
##         : $script_path             => Path to variant_effect_predictor script
##         : $assembly                => Assembly version to use
##         : $cache_directory         => VEP chache directory
##         : $buffer_size             => Sets the internal buffer size, corresponding to the number of variants that are read in to memory simultaneously
##         : $infile_format           => Input file format - one of "ensembl", "vcf", "hgvs", "id"
##         : $outfile_format          => Output file format

    my ($arg_href) = @_;

    ## Default(s)
    my $infile_format;
    my $outfile_format;
    my $fork;

    ## Flatten argument(s)
    my $plugins_ref;
    my $regions_ref;
    my $vep_features_ref;
    my $outfile_path;
    my $infile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $reference_path;
    my $script_path;
    my $assembly;
    my $cache_directory;
    my $buffer_size;

    my $tmpl = {
        plugins_ref =>
          { default => [], strict_type => 1, store => \$plugins_ref },
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        vep_features_ref =>
          { default => [], strict_type => 1, store => \$vep_features_ref },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        infile_path     => { strict_type => 1, store => \$infile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE      => { store       => \$FILEHANDLE },
        reference_path  => { strict_type => 1, store => \$reference_path },
        script_path     => { strict_type => 1, store => \$script_path },
        assembly        => { strict_type => 1, store => \$assembly },
        cache_directory => { strict_type => 1, store => \$cache_directory },
        buffer_size => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$buffer_size
        },
        infile_format => {
            default     => q{vcf},
            allow       => [qw{ ensembl vcf hgvs id }],
            strict_type => 1,
            store       => \$infile_format
        },
        outfile_format => {
            default     => q{vcf},
            allow       => [qw{ vcf tab json }],
            strict_type => 1,
            store       => \$outfile_format
        },
        fork => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$fork
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Vep
    # Stores commands depending on input parameters
    my @commands = qw{ perl };

    if ($script_path) {

        push @commands, $script_path;
    }
    else {

        push @commands, q{vep};
    }

    ## Options
    if ($fork) {

        push @commands, q{--fork} . $SPACE . $fork;
    }
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
            stdoutfile_path        => $stdoutfile_path,
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

sub variant_effect_predictor_install {

## variant_effect_predictor_install

## Function : Perl wrapper for vep INSTALL script. Based on version 90.
## Returns  : @commands

## Arguments: $plugins_ref, $species_ref, $auto, $cache_directory, $assembly, $stdoutfile_path, $stderrfile_path, stderrfile_path_append, $FILEHANDLE
##          : $plugins_ref            => Vep plugins {REF}
##          : $species_ref            => Comma-separated list of species to install when using --AUTO {REF}
##          : $auto                   => Run installer without user prompts. Use "a" (API + Faidx/htslib),"l" (Faidx/htslib only), "c" (cache), "f" (FASTA), "p" (plugins) to specify parts to install.
##          : $cache_directory        => Set destination directory for cache files
##          : $assembly               => Assembly name to use if more than one during --AUTO
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $plugins_ref;
    my $species_ref;
    my $auto;
    my $cache_directory;
    my $assembly;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        plugins_ref =>
          { default => [], strict_type => 1, store => \$plugins_ref },
        species_ref => {
            default     => [qw{ homo_sapiens }],
            strict_type => 1,
            store       => \$species_ref
        },
        auto            => { strict_type => 1, store => \$auto },
        cache_directory => { strict_type => 1, store => \$cache_directory },
        assembly        => { strict_type => 1, store => \$assembly },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store => \$FILEHANDLE },

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
            stdoutfile_path        => $stdoutfile_path,
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
