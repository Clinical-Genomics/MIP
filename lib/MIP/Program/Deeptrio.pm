package MIP::Program::Deeptrio;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $EQUALS $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ deeptrio };
}

sub deeptrio {

## Function : Perl wrapper for generic commands module
## Returns  : @commands
## Arguments: $bedfile                => Bed file containing the list of regions we want to process
##          : $filehandle             => Filehandle to write to
##          : $model_type             => Type of model to use for variant calling. Allowed values WES, WGS, or PACBIO
##          : $num_shards             => Number of files the input is split into for the make examples step
##          : $output_gvcf_child      => Path to output gvcf file for the child in a duo or trio
##          : $output_gvcf_parent1    => Path to output gvcf for a parent in a duo or a trio
##          : $output_gvcf_parent2    => Path to output gvcf for a parent in trio
##          : $output_vcf_child       => Path to output vcf file for the child in a duo or trio
##          : $output_vcf_parent1     => Path to output vcf for a parent in a duo or trio
##          : $output_vcf_parent2     => Path to output vcf for a parent in trio
##          : $reads_child            => Path to input bam for the child in a duo or trio
##          : $reads_parent1          => Path to input bam for a parent in a duo or trio
##          : $reads_parent2          => Path to input bam for a parent in trio
##          : $referencefile_path     => Path to genome reference
##          : $sample_name_child      => Sample name for the child in a duo or trio
##          : $sample_name_parent1    => Sample name for a parent in a duo or trio
##          : $sample_name_parent2    => Sample name for a parent in trio
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bedfile;
    my $filehandle;
    my $model_type;
    my $num_shards;
    my $output_gvcf_child;
    my $output_gvcf_parent1;
    my $output_gvcf_parent2;
    my $output_vcf_child;
    my $output_vcf_parent1;
    my $output_vcf_parent2;
    my $reads_child;
    my $reads_parent1;
    my $reads_parent2;
    my $referencefile_path;
    my $sample_name_child;
    my $sample_name_parent1;
    my $sample_name_parent2;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    my $tmpl = {
        bedfile => {
            defined     => 1,
            store       => \$bedfile,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        model_type => {
            allow       => [qw{ WES WGS PACBIO }],
            defined     => 1,
            required    => 1,
            store       => \$model_type,
            strict_type => 1,
        },
        num_shards => {
            allow       => qr/ \A \d+ \z /xms,
            required    => 1,
            store       => \$num_shards,
            strict_type => 1,
        },
        output_gvcf_child => {
            defined     => 1,
            required    => 1,
            store       => \$output_gvcf_child,
            strict_type => 1,
        },
        output_gvcf_parent1 => {
            defined     => 1,
            required    => 1,
            store       => \$output_gvcf_parent1,
            strict_type => 1,
        },
        output_gvcf_parent2 => {
            store       => \$output_gvcf_parent2,
            strict_type => 1,
        },
        output_vcf_child => {
            defined     => 1,
            required    => 1,
            store       => \$output_vcf_child,
            strict_type => 1,
        },
        output_vcf_parent1 => {
            defined     => 1,
            required    => 1,
            store       => \$output_vcf_parent1,
            strict_type => 1,
        },
        output_vcf_parent2 => {
            store       => \$output_vcf_parent2,
            strict_type => 1,
        },
        reads_child => {
            defined     => 1,
            required    => 1,
            store       => \$reads_child,
            strict_type => 1,
        },
        reads_parent1 => {
            defined     => 1,
            required    => 1,
            store       => \$reads_parent1,
            strict_type => 1,
        },
        reads_parent2 => {
            store       => \$reads_parent2,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        sample_name_child => {
            defined     => 1,
            required    => 1,
            store       => \$sample_name_child,
            strict_type => 1,
        },
        sample_name_parent1 => {
            defined     => 1,
            required    => 1,
            store       => \$sample_name_parent1,
            strict_type => 1,
        },
        sample_name_parent2 => {
            store       => \$sample_name_parent2,
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
        stdinfile_path => {
            store       => \$stdinfile_path,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => q{run_deeptrio}, } ), );

    push @commands, q{--model_type} . $SPACE . $model_type;
    push @commands, q{--ref} . $SPACE . $referencefile_path;
    push @commands, q{--num_shards} . $SPACE . $num_shards;

    ## Child
    push @commands, q{--sample_name_child} . $SPACE . $sample_name_child;
    push @commands, q{--reads_child} . $SPACE . $reads_child;
    push @commands, q{--output_gvcf_child} . $SPACE . $output_gvcf_child;
    push @commands, q{--output_vcf_child} . $SPACE . $output_vcf_child;

    ## Parent1
    push @commands, q{--sample_name_parent1} . $SPACE . $sample_name_parent1;
    push @commands, q{--reads_parent1} . $SPACE . $reads_parent1;
    push @commands, q{--output_gvcf_parent1} . $SPACE . $output_gvcf_parent1;
    push @commands, q{--output_vcf_parent1} . $SPACE . $output_vcf_parent1;

    if ($sample_name_parent2) {
        push @commands, q{--sample_name_parent2} . $SPACE . $sample_name_parent2;
    }
    if ($reads_parent2) {
        push @commands, q{--reads_parent2} . $SPACE . $reads_parent2;
    }
    if ($output_gvcf_parent2) {
        push @commands, q{--output_gvcf_parent2} . $SPACE . $output_gvcf_parent2;
    }
    if ($output_vcf_parent2) {
        push @commands, q{--output_vcf_parent2} . $SPACE . $output_vcf_parent2;
    }
    if ($bedfile) {
        push @commands, q{--regions} . $SPACE . $bedfile;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
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
