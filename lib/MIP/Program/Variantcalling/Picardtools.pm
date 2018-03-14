package MIP::Program::Variantcalling::Picardtools;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Language::Java qw{java_core};
use MIP::Program::Base::Picardtools qw{ picardtools_base};
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ picardtools_genotypeconcordance picardtools_sortvcf sort_vcf };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub picardtools_sortvcf {

## Function : Perl wrapper for writing picardtools sortvcf recipe to $FILEHANDLE. Based on picardtools 2.14.1.
## Returns  : @commands
## Arguments: $create_index         => Create index
##          : $FILEHANDLE           => Sbatch filehandle to write to
##          : $java_jar             => Java jar
##          : $java_use_large_pages => Use java large pages
##          : $infile_paths_ref     => Infile paths {REF}
##          : $memory_allocation    => Memory allocation for java
##          : $outfile_path         => Outfile path
##          : $referencefile_path   => Genome reference file
##          : $sequence_dictionary  => Sequence dictionary
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $referencefile_path;
    my $sequence_dictionary;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $create_index;
    my $java_use_large_pages;

    my $tmpl = {
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
        FILEHANDLE       => { store => \$FILEHANDLE },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        java_jar             => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        outfile_path      => {
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
        temp_directory  => { strict_type => 1, store => \$temp_directory },
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

## Function : Perl wrapper for writing picardtools genotypeconcordance recipe to $FILEHANDLE. Based on picardtools 2.14.1.
## Returns  : @commands
## Arguments: $call_sample          => Name of the call sample within the call VCF
##          : $create_index         => Create index
##          : $FILEHANDLE           => Sbatch filehandle to write to
##          : $infile_path          => Infile paths
##          : $intervals_ref        => One or more genomic intervals over which to operate {REF}
##          : $java_jar             => Java jar
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memory allocation for java
##          : $min_depth            => Genotypes below this depth will have genotypes classified as LowDp
##          : $min_genotype_quality => Genotypes below this genotype quality will have genotypes classified as LowGq
##          : $outfile_prefix_path  => Outfile path
##          : $referencefile_path   => Genome reference file
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $truth_file_path      => VCF containing the truth sample
##          : $truth_sample         => Name of the truth sample within the truth VCF

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $intervals_ref;
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
    my $create_index;
    my $java_use_large_pages;
    my $min_depth;
    my $min_genotype_quality;

    my $tmpl = {
        call_sample => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$call_sample
        },
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
        FILEHANDLE  => { store => \$FILEHANDLE },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        java_jar             => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        min_depth         => {
            default     => 0,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$min_depth
        },
        min_genotype_quality => {
            default     => 0,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$min_genotype_quality
        },
        outfile_prefix_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_prefix_path
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        truth_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$truth_file_path
        },
        truth_sample => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$truth_sample
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        temp_directory  => { strict_type => 1, store => \$temp_directory },
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

sub sort_vcf {

## Function : Writes sbatch code to supplied filehandle to sort variants in vcf format.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $FILEHANDLE            => SBATCH script FILEHANDLE to print to
##          : $infile_paths_ref      => Infiles to sort {REF}
##          : $outfile               => The sorted outfile
##          : $sequence_dict_file    => Human reference sequence dict file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $outfile;
    my $sequence_dict_file;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE, },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        outfile => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile,
        },
        sequence_dict_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sequence_dict_file
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {$FILEHANDLE} q{## Picard SortVcf};

    ## Writes java core commands to filehandle.
    picardtools_sortvcf(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => \@{$infile_paths_ref},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            java_jar => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            memory_allocation => q{Xmx2g},
            outfile_path      => $outfile,
            referencefile_path =>
              $active_parameter_href->{human_genome_reference},
            sequence_dictionary => $sequence_dict_file,
            temp_directory      => $active_parameter_href->{temp_directory},
        }
    );
    return;
}

1;
