package MIP::Program::Variantcalling::Picardtools;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };
use MIP::Language::Java qw{java_core};
use MIP::Program::Base::Picardtools qw{ picardtools_base};

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ picardtools_sortvcf picardtools_genotypeconcordance };
}

## Constants
Readonly my $SPACE => q{ };

sub picardtools_sortvcf {

## Function : Perl wrapper for writing picardtools sortvcf recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $infile_paths_ref     => Infile paths {REF}
##          : $outfile_path         => Outfile path
##          : $referencefile_path   => Genome reference file
##          : $sequence_dictionary  => Sequence dictionary
##          : $FILEHANDLE           => Sbatch filehandle to write to
##          : $stderrfile_path      => Stderrfile path
##          : $memory_allocation    => Memory allocation for java
##          : $temp_directory       => Redirect tmp files to java temp
##          : $java_use_large_pages => Use java large pages
##          : $java_jar             => Java jar
##          : $create_index         => Create index

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $sequence_dictionary;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $memory_allocation;
    my $temp_directory;
    my $java_jar;

    ## Default(s)
    my $java_use_large_pages;
    my $create_index;

    my $tmpl = {
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        sequence_dictionary => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sequence_dictionary
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        temp_directory    => { strict_type => 1, store => \$temp_directory },
        java_jar          => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
    if ($java_jar) {

        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }

    ## Picardtools mergesamfiles
    push @commands, q{SortVcf};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    push @commands, q{SEQUENCE_DICTIONARY=} . $sequence_dictionary;

    ## Infile
    push @commands, q{INPUT=} . join $SPACE . q{INPUT=}, @{$infile_paths_ref};

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

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

sub picardtools_genotypeconcordance {

## Function : Perl wrapper for writing picardtools genotypeconcordance recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $intervals_ref        => One or more genomic intervals over which to operate {REF}
##          : $infile_path          => Infile paths
##          : $truth_file_path      => VCF containing the truth sample
##          : $outfile_prefix_path  => Outfile path
##          : $truth_sample         => Name of the truth sample within the truth VCF
##          : $call_sample          => Name of the call sample within the call VCF
##          : $referencefile_path   => Genome reference file
##          : $FILEHANDLE           => Sbatch filehandle to write to
##          : $stderrfile_path      => Stderrfile path
##          : $memory_allocation    => Memory allocation for java
##          : $temp_directory       => Redirect tmp files to java temp
##          : $java_use_large_pages => Use java large pages
##          : $java_jar             => Java jar
##          : $create_index         => Create index
##          : $min_genotype_quality => Genotypes below this genotype quality will have genotypes classified as LowGq
##          : $min_depth            => Genotypes below this depth will have genotypes classified as LowDp

    my ($arg_href) = @_;

    ## Default(s)
    my $min_genotype_quality;
    my $min_depth;

    ## Flatten argument(s)
    my $intervals_ref;
    my $infile_path;
    my $truth_file_path;
    my $outfile_prefix_path;
    my $truth_sample;
    my $call_sample;
    my $referencefile_path;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $memory_allocation;
    my $temp_directory;
    my $java_jar;

    ## Default(s)
    my $java_use_large_pages;
    my $create_index;

    my $tmpl = {
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        truth_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$truth_file_path
        },
        outfile_prefix_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_prefix_path
        },
        truth_sample => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$truth_sample
        },
        call_sample => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$call_sample
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        stderrfile_path   => { strict_type => 1, store => \$stderrfile_path },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        temp_directory    => { strict_type => 1, store => \$temp_directory },
        java_jar          => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
        min_genotype_quality => {
            default     => 0,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$min_genotype_quality
        },
        min_depth => {
            default     => 0,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$min_depth
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
    if ($java_jar) {

        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }

    ## Picardtools mergesamfiles
    push @commands, q{GenotypeConcordance};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    ##Options
    if ($min_genotype_quality) {

        push @commands, q{MIN_GQ=} . $min_genotype_quality;
    }
    if ($min_depth) {

        push @commands, q{MIN_DP=} . $min_depth;
    }
    if ( @{$intervals_ref} ) {

        push @commands, q{INTERVALS=} . join $SPACE . q{INTERVALS=},
          @{$intervals_ref};
    }

    ## Infile
    push @commands, q{CALL_VCF=} . $infile_path;

    push @commands, q{TRUTH_VCF=} . $truth_file_path;

    ## Output
    push @commands, q{OUTPUT=} . $outfile_prefix_path;

    push @commands, q{TRUTH_SAMPLE=} . $truth_sample;

    push @commands, q{CALL_SAMPLE=} . $call_sample;

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
