package MIP::Program::Delly;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ delly_call delly_filter delly_merge };
}

## Constants
Readonly my $MIN_MAP_QUAL       => 20;
Readonly my $INSERT_SIZE_CUTOFF => 15;

sub delly_call {

## Function : Perl wrapper for writing Delly call recipe to $filehandle or return commands array. Based on Delly 0.7.8.
## Returns  : @commands
## Arguments: $exclude_file_path      => File with regions to exclude
##          : $filehandle             => Filehandle to write to
##          : $genotypefile_path      => Input VCF/BCF file for re-genotyping
##          : $infile_path            => Infile path
##          : $mad_cutoff             => Insert size cutoff, median+s*MAD (deletions only)
##          : $mapping_qual           => Minimum paired-end mapping quality
##          : $small_indel            => Small indel calling
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Reference sequence file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $sv_type                => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exclude_file_path;
    my $filehandle;
    my $genotypefile_path;
    my $infile_path;
    my $mad_cutoff;
    my $mapping_qual;
    my $small_indel;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $sv_type;

    my $tmpl = {
        exclude_file_path => { store => \$exclude_file_path, strict_type => 1, },
        filehandle        => {
            store => \$filehandle,
        },
        genotypefile_path => { store => \$genotypefile_path, strict_type => 1, },
        infile_path       => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        mad_cutoff => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            default     => $INSERT_SIZE_CUTOFF,
            store       => \$mad_cutoff,
            strict_type => 1,
        },
        mapping_qual => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            default     => $MIN_MAP_QUAL,
            store       => \$mapping_qual,
            strict_type => 1,
        },
        small_indel => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$small_indel,
            strict_type => 1,
        },
        outfile_path       => { store => \$outfile_path, strict_type => 1, },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
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
        sv_type => {
            allow       => [qw{ DEL DUP INV INS BND }],
            store       => \$sv_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ delly call };

    if ($sv_type) {

        push @commands, q{--type} . $SPACE . $sv_type;
    }

    push @commands, q{--map-qual} . $SPACE . $mapping_qual;

    push @commands, q{--mad-cutoff} . $SPACE . $mad_cutoff;

    if ($small_indel) {

        push @commands, q{--i};
    }

    if ($exclude_file_path) {

        push @commands, q{--exclude} . $SPACE . $exclude_file_path;
    }

    if ($referencefile_path) {

        push @commands, q{--genome} . $SPACE . $referencefile_path;
    }

    if ($genotypefile_path) {

        push @commands, q{--vcffile} . $SPACE . $genotypefile_path;
    }

    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    push @commands, $infile_path;

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

sub delly_merge {

## Function : Perl wrapper for writing Delly merge recipe to $filehandle or return commands array. Based on Delly 0.7.8.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_paths_ref       => Infile paths {REF}
##          : $max_size               => Max. SV size
##          : $min_size               => Min. SV size
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $sv_type                => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $max_size;
    my $min_size;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $sv_type;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref,
        },
        max_size => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            strict_type => 1,
            store       => \$max_size,
        },
        min_size => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            strict_type => 1,
            store       => \$min_size,
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
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
        sv_type => {
            allow       => [qw{ DEL DUP INV INS BND }],
            strict_type => 1,
            store       => \$sv_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ delly merge };

    if ($sv_type) {

        push @commands, q{--type} . $SPACE . $sv_type;
    }

    if ($min_size) {

        push @commands, q{--minsize} . $SPACE . $min_size;
    }

    if ($max_size) {

        push @commands, q{--maxsize} . $SPACE . $max_size;
    }

    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

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

sub delly_filter {

## Function : Perl wrapper for writing Delly filter recipe to $filehandle or return commands array. Based on Delly 0.7.6.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $filter_mode            => Filter mode
##          : $infile_path            => Infile pat
##          : $max_size               => Max. SV size
##          : $min_size               => Min. SV size
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $sv_type                => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $filter_mode;
    my $infile_path;
    my $max_size;
    my $min_size;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $sv_type;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        filter_mode => {
            required    => 1,
            defined     => 1,
            allow       => [qw{ somatic germline }],
            strict_type => 1,
            store       => \$filter_mode,
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        max_size => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            strict_type => 1,
            store       => \$max_size,
        },
        min_size => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            strict_type => 1,
            store       => \$min_size,
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
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
        sv_type => {
            required    => 1,
            defined     => 1,
            allow       => [qw{ DEL DUP INV INS BND }],
            strict_type => 1,
            store       => \$sv_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ delly filter };

    if ($sv_type) {

        push @commands, q{--type} . $SPACE . $sv_type;
    }

    if ($filter_mode) {

        push @commands, q{--filter} . $SPACE . $filter_mode;
    }

    if ($min_size) {

        push @commands, q{--minsize} . $SPACE . $min_size;
    }

    if ($max_size) {

        push @commands, q{--maxsize} . $SPACE . $max_size;
    }

    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    push @commands, $infile_path;

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
