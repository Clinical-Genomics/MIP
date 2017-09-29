package MIP::Program::Variantcalling::Cnvnator;

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
use File::Spec::Functions qw{ catdir };

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
    our @EXPORT_OK =
      qw{ cnvnator_read_extraction cnvnator_histogram cnvnator_statistics cnvnator_partition cnvnator_calling cnvnator_convert_to_vcf };

}

## Constants
Readonly my $SPACE => q{ };

sub cnvnator_read_extraction {

## cnvnator_read_extraction

## Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $regions_ref, $infile_paths_ref, $outfile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $unique
##          : $regions_ref      => The regions to process {REF}
##          : $infile_paths_ref => Infile paths {REF}
##          : $outfile_path     => Outfile path
##          : $stderrfile_path  => Stderrfile path
##          : $stdoutfile_path  => Stdoutfile path
##          : $FILEHANDLE       => Filehandle to write to
##          : $unique           => Ensure correct q0 field for CNV calls

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $unique;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
        unique     => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$unique
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options
    if ( @{$regions_ref} ) {    #Limit output to regions

        push @commands, q{-chrom}, join $SPACE, @{$regions_ref};
    }

    if ($unique) {

        push @commands, q{-unique};
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{-root} . $SPACE . $outfile_path;
    }

    ## Infile
    push @commands, q{-tree}, join $SPACE, @{$infile_paths_ref};

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
            stdoutfile_path => $stdoutfile_path,
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

sub cnvnator_histogram {

## cnvnator_histogram

## Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $regions_ref, $infile_path, $referencedirectory_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $cnv_bin_size
##          : $regions_ref             => The regions to process {REF}
##          : $infile_path             => Infile paths
##          : $referencedirectory_path => Reference sequence file
##          : $stderrfile_path         => Stderrfile path
##          : $stdoutfile_path         => Stdoutfile path
##          : $FILEHANDLE              => Filehandle to write to
##          : $cnv_bin_size            => Copy number variant bin size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $referencedirectory_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $cnv_bin_size;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        referencedirectory_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencedirectory_path
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE   => { store => \$FILEHANDLE },
        cnv_bin_size => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$cnv_bin_size
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options: mimit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{-chrom}, join $SPACE, @{$regions_ref};
    }

    if ($referencedirectory_path) {

        push @commands, q{-d} . $SPACE . $referencedirectory_path;
    }

    if ($cnv_bin_size) {

        push @commands, q{-his} . $SPACE . $cnv_bin_size;
    }

    ## Infile
    push @commands, q{-root}, $infile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
            stdoutfile_path => $stdoutfile_path,
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

sub cnvnator_statistics {

## cnvnator_statistics

## Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $regions_ref, $infile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $cnv_bin_size
##          : $regions_ref     => The regions to process {REF}
##          : $infile_path     => Infile paths
##          : $stderrfile_path => Stderrfile path
##          : $stdoutfile_path => Stdoutfile path
##          : $FILEHANDLE      => Filehandle to write to
##          : $cnv_bin_size    => Copy number variant bin size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $cnv_bin_size;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE   => { store => \$FILEHANDLE },
        cnv_bin_size => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$cnv_bin_size
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options: limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{-chrom}, join $SPACE, @{$regions_ref};
    }

    if ($cnv_bin_size) {

        push @commands, q{-stat} . $SPACE . $cnv_bin_size;
    }

    ## Infile
    push @commands, q{-root}, $infile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
            stdoutfile_path => $stdoutfile_path,
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

sub cnvnator_partition {

## cnvnator_partition

## Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $regions_ref, $infile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $cnv_bin_size
##          : $regions_ref     => The regions to process {REF}
##          : $infile_path     => Infile paths
##          : $stderrfile_path => Stderrfile path
##          : $stdoutfile_path => Stdoutfile path
##          : $FILEHANDLE      => Filehandle to write to
##           : $cnv_bin_size    => Copy number variant bin size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $cnv_bin_size;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE   => { store => \$FILEHANDLE },
        cnv_bin_size => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$cnv_bin_size
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options: limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{-chrom}, join $SPACE, @{$regions_ref};
    }

    if ($cnv_bin_size) {

        push @commands, q{-partition} . $SPACE . $cnv_bin_size;
    }

    ## Infile
    push @commands, q{-root}, $infile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
            stdoutfile_path => $stdoutfile_path,
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

sub cnvnator_calling {

## cnvnator_calling

## Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $regions_ref, $infile_path, $outfile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $cnv_bin_size
##          : $regions_ref             => The regions to process {REF}
##          : $infile_path             => Infile paths
##          : $outfile_path            => Outfile path
##          : $stderrfile_path         => Stderrfile path
##          : $stdoutfile_path         => Stdoutfile path
##          : $FILEHANDLE              => Filehandle to write to
##          : $cnv_bin_size            => Copy number variant bin size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $cnv_bin_size;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE   => { store => \$FILEHANDLE },
        cnv_bin_size => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$cnv_bin_size
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options: limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{-chrom}, join $SPACE, @{$regions_ref};
    }

    if ($cnv_bin_size) {

        push @commands, q{-call} . $SPACE . $cnv_bin_size;
    }

    ## Infile
    push @commands, q{-root}, $infile_path;

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{>} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
            stdoutfile_path => $stdoutfile_path,
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

sub cnvnator_convert_to_vcf {

## cnvnator_convert_to_vcf

## Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $infile_path, $referencedirectory_path, $outfile_path, $stderrfile_path, $FILEHANDLE
##          : $infile_path             => Infile paths
##          : $referencedirectory_path => Reference sequence file
##          : $outfile_path            => Stdoutfile path
##          : $stderrfile_path         => Stderrfile path
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $referencedirectory_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        referencedirectory_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencedirectory_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator2VCF.pl };

    ## Infile
    push @commands, $infile_path;

    ## Options
    if ($referencedirectory_path) {

        push @commands, $referencedirectory_path;
    }

    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
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
