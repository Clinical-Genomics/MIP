package MIP::Program::Rtg;

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
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ rtg_format rtg_vcfeval };
}

sub rtg_format {

## Function : Perl wrapper for rtg tools 3.9.1.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $input_format           => Format of input
##          : $reference_genome_path  => Human reference genome file path
##          : $sdf_output_directory   => Directory name of output SDF
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $reference_genome_path;
    my $sdf_output_directory;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $input_format;

    my $tmpl = {
        input_format => {
            allow       => [qw{ fasta fastq fastq-interleaved sam-se sam-pe }],
            default     => q{fasta},
            store       => \$input_format,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        reference_genome_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_path,
            strict_type => 1,
        },
        sdf_output_directory => {
            defined     => 1,
            required    => 1,
            store       => \$sdf_output_directory,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ rtg format };

    push @commands, q{--format=} . $input_format;

    push @commands, q{--output=} . $sdf_output_directory;

    push @commands, $reference_genome_path;

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

sub rtg_vcfeval {

## Function : Perl wrapper for rtg tools 3.9.1.
## Returns  : @commands
## Arguments: $all_record             => Use all records regardless of FILTER status
##          : $bed_regionsfile_path   => If set, only read VCF records that overlap the ranges contained in the specified BED file
##          : $baselinefile_path      => VCF file containing baseline variants
##          : $callfile_path          => VCF file containing called variants
##          : $eval_region_file_path  => Evaluate within regions contained in the supplied BED file, allowing transborder matches
##          : $FILEHANDLE             => Filehandle to write to
##          : $outputdirectory_path   => Directory for output
##          : $output_mode            => Output reporting mode
##          : $sample_id              => Sample ID
##          : $sdf_template_file_path => SDF (SDF=Rtg specif format) of the reference genome the variants are called against
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $thread_number          => Number of threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $baselinefile_path;
    my $bed_regionsfile_path;
    my $callfile_path;
    my $eval_region_file_path;
    my $FILEHANDLE;
    my $outputdirectory_path;
    my $sample_id;
    my $sdf_template_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $thread_number;

    ## Default(s)
    my $all_record;
    my $output_mode;

    my $tmpl = {
        all_record => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$all_record,
            strict_type => 1,
        },
        bed_regionsfile_path => {
            store       => \$bed_regionsfile_path,
            strict_type => 1,
        },
        baselinefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$baselinefile_path,
            strict_type => 1,
        },
        callfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$callfile_path,
            strict_type => 1,
        },
        eval_region_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$eval_region_file_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        outputdirectory_path => {
            defined     => 1,
            required    => 1,
            store       => \$outputdirectory_path,
            strict_type => 1,
        },
        output_mode => {
            allow       => [qw{ split annotate combine ga4gh roc-only }],
            default     => q{split},
            store       => \$output_mode,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
        sdf_template_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$sdf_template_file_path,
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
        thread_number => {
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ rtg vcfeval };

    if ($all_record) {

        push @commands, q{--all-records};
    }

    if ($thread_number) {

        push @commands, q{--threads=} . $thread_number;
    }

    push @commands, q{--baseline=} . $baselinefile_path;

    if ($bed_regionsfile_path) {

        push @commands, q{--bed-regions=} . $bed_regionsfile_path;
    }
    push @commands, q{--calls=} . $callfile_path;

    push @commands, q{--evaluation-regions=} . $eval_region_file_path;

    push @commands, q{--template=} . $sdf_template_file_path;

    if ($sample_id) {

        push @commands, q{--sample=} . $sample_id;
    }

    push @commands, q{--output-mode=} . $output_mode;

    push @commands, q(--output=) . $outputdirectory_path;

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
