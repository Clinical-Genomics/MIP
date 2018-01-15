package MIP::Program::Qc::Contest;

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
use MIP::Program::Base::Gatk qw{ gatk_base };
use MIP::Language::Java qw{ java_core };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gatk_contest };
}

## Constants
Readonly my $SPACE => q{ };

sub gatk_contest {

## Function : Perl wrapper for writing GATK ContEst recipe to $FILEHANDLE. Based on GATK 3.8.0.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $referencefile_path     => Reference sequence file
##          : $outfile_path           => Outfile path
##          : $popvcffilepath         => Population allele frequency vcf file path
##          : $infile_eval            => Input bam file for evaluation
##          : $infile_genotype        => Input bam file for genotyping on the fly
##          : $min_genotype_ratio     => The ratio of alt to other bases
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $downsample_to_coverage;
    my $FILEHANDLE;
    my $infile_eval;
    my $infile_genotype;
    my $intervals_ref;
    my $java_use_large_pages;
    my $java_jar;
    my $memory_allocation;
    my $min_genotype_ratio;
    my $outfile_path;
    my $pedigree;
    my $popvcffilepath;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temp_directory;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;

    ## Default(s)
    ## None

    my $tmpl = {
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        gatk_disable_auto_index_and_file_lock => {
            default     => 0,
            allow       => [ 0, 1 ],
            store       => \$gatk_disable_auto_index_and_file_lock,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            store       => \$intervals_ref,
            strict_type => 1,
        },
        infile_eval => {
            defined     => 1,
            required    => 1,
            store       => \$infile_eval,
            strict_type => 1,
        },
        infile_genotype => {
            defined     => 1,
            required    => 1,
            store       => \$infile_genotype,
            strict_type => 1,
        },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        java_jar => {
            store       => \$java_jar,
            strict_type => 1,
        },
        logging_level => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$logging_level,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        min_genotype_ratio => {
            ## Exactly 2 decimal points after 0 or 1
            allow       => qr/ ^0.\d{1,2}$ | ^1$ /xsm,
            defined     => 1,
            default     => 0.01,
            required    => 0,
            store       => \$min_genotype_ratio,
            strict_type => 1,
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree => {
            store       => \$pedigree,
            strict_type => 1,
        },
        pedigree_validation_type => {
            allow       => [qw{ STRICT SILENT }],
            default     => q{SILENT},
            store       => \$pedigree_validation_type,
            strict_type => 1,
        },
        popvcffilepath => {
            defined     => 1,
            required    => 1,
            store       => \$popvcffilepath,
            strict_type => 1,
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
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
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands;

    if ($java_jar) {    # Write java core commands to filehandle.
        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            commands_ref             => \@commands,
            analysis_type            => q{ContEst},
            logging_level            => $logging_level,
            intervals_ref            => $intervals_ref,
            referencefile_path       => $referencefile_path,
            pedigree                 => $pedigree,
            pedigree_validation_type => $pedigree_validation_type,
            downsample_to_coverage   => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
        }
    );

    ############################################
    ## ADD COMMAND SPECIFIC FLAGS AND OPTIONS ##
    ############################################

    if ($min_genotype_ratio) {
        push @commands, q{--min_genotype_ratio} . $SPACE . $min_genotype_ratio;
    }

    push @commands, q{--input_file:eval} . $SPACE .  $infile_eval;

    push @commands, q{--input_file:genotype} . $SPACE .  $infile_genotype;

    push @commands, q{--popfile} . $SPACE . $popvcffilepath;

    push @commands, q{--out} . $SPACE . $outfile_path;

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
