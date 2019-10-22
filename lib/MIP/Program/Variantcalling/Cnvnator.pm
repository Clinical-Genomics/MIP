package MIP::Program::Variantcalling::Cnvnator;

use 5.026;
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
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ cnvnator_read_extraction cnvnator_histogram cnvnator_statistics cnvnator_partition cnvnator_calling cnvnator_convert_to_vcf };

}

## Constants
Readonly my $SPACE => q{ };

sub cnvnator_read_extraction {

## Function : Perl wrapper for writing cnvnator recipe to $filehandle or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $filehandle       => Filehandle to write to
##          : $infile_paths_ref => Infile paths {REF}
##          : $outfile_path      => outfile_path
##          : $regions_ref      => The regions to process {REF}
##          : $stderrfile_path  => Stderrfile path
##          : $stdoutfile_path  => Stdoutfile path
##          : $unique           => Ensure correct q0 field for CNV calls

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $outfile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $unique;

    my $tmpl = {
        filehandle       => { store => \$filehandle, },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        unique          => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$unique,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options
    # Limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{-chrom} . $SPACE . join $SPACE, @{$regions_ref};
    }

    if ($unique) {

        push @commands, q{-unique};
    }

    if ($outfile_path) {

        # Specify output filename
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );

    return @commands;

}

sub cnvnator_histogram {

## Function : Perl wrapper for writing cnvnator recipe to $filehandle or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $cnv_bin_size            => Copy number variant bin size
##          : $filehandle              => Filehandle to write to
##          : $infile_path             => Infile paths
##          : $referencedirectory_path => Reference sequence file
##          : $regions_ref             => The regions to process {REF}
##          : $stderrfile_path         => Stderrfile path
##          : $stdoutfile_path         => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cnv_bin_size;
    my $filehandle;
    my $infile_path;
    my $referencedirectory_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stdoutfile_path;

    my $tmpl = {
        cnv_bin_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$cnv_bin_size,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        referencedirectory_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencedirectory_path,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options:
    # limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{-chrom} . $SPACE . join $SPACE, @{$regions_ref};
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );

    return @commands;

}

sub cnvnator_statistics {

## Function : Perl wrapper for writing cnvnator recipe to $filehandle or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $cnv_bin_size    => Copy number variant bin size
##          : $filehandle      => Filehandle to write to
##          : $infile_path     => Infile paths
##          : $regions_ref     => The regions to process {REF}
##          : $stderrfile_path => Stderrfile path
##          : $stdoutfile_path => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cnv_bin_size;
    my $filehandle;
    my $infile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle   => { store => \$filehandle, },
        cnv_bin_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$cnv_bin_size,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options:
    # Limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{-chrom} . $SPACE . join $SPACE, @{$regions_ref};
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );

    return @commands;

}

sub cnvnator_partition {

## Function : Perl wrapper for writing cnvnator recipe to $filehandle or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $cnv_bin_size    => Copy number variant bin size
##          : $filehandle      => Filehandle to write to
##          : $infile_path     => Infile paths
##          : $regions_ref     => The regions to process {REF}
##          : $stderrfile_path => Stderrfile path
##          : $stdoutfile_path => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cnv_bin_size;
    my $filehandle;
    my $infile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle   => { store => \$filehandle, },
        cnv_bin_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$cnv_bin_size,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options:
    # Limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{-chrom} . $SPACE . join $SPACE, @{$regions_ref};
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );

    return @commands;

}

sub cnvnator_calling {

## Function : Perl wrapper for writing cnvnator recipe to $filehandle or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $cnv_bin_size    => Copy number variant bin size
##          : $filehandle      => Filehandle to write to
##          : $infile_path     => Infile paths
##          : $regions_ref     => The regions to process {REF}
##          : $stderrfile_path => Stderrfile path
##          : $stdoutfile_path => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cnv_bin_size;
    my $filehandle;
    my $infile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stdoutfile_path;

    my $tmpl = {
        cnv_bin_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$cnv_bin_size,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## cnvnator
    my @commands = qw{ cnvnator };

    ## Options:
    # Limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{-chrom} . $SPACE . join $SPACE, @{$regions_ref};
    }

    if ($cnv_bin_size) {

        push @commands, q{-call} . $SPACE . $cnv_bin_size;
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );

    return @commands;

}

sub cnvnator_convert_to_vcf {

## Function : Perl wrapper for writing cnvnator recipe to $filehandle or return commands array. Based on cnvnator 0.3.3.
## Returns  : "@commands"
## Arguments: $filehandle              => Filehandle to write to
##          : $infile_path             => Infile paths
##          : $referencedirectory_path => Reference sequence file
##          : $stdoutfile_path         => Stdoutfile path
##          : $stderrfile_path         => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $referencedirectory_path;
    my $stdoutfile_path;
    my $stderrfile_path;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        referencedirectory_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencedirectory_path,
            strict_type => 1,
        },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );

    return @commands;

}

1;
