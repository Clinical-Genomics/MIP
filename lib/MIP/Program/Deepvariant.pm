package MIP::Program::Deepvariant;

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
use MIP::Constants qw{ $SPACE $EQUALS };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ deepvariant };
}

sub deepvariant {

## Function : Perl wrapper for generic commands module
## Returns  : @commands
## Arguments: $bedfile                  => Bed file containing the list of regions we want to process
##          : $filehandle               => Filehandle to write to
##          : $infile_path              => Aligned, sorted, indexed bam file containing the reads we want to call
##          : $intermediate_results_dir => Store intermediate results here
##          : $model_type               => Type of model to use for variant calling. Allowed values WES, WGS, or PACBIO
##          : $num_shards               => Number of files the input is split into for the make examples step
##          : $outfile_path             => Path to the output gvcf file
##          : $outfile_path_vcf         => Path to the output vcf file
##          : $referencefile_path       => Path to genome reference
##          : $stderrfile_path          => Stderrfile path
##          : $stderrfile_path_append   => Append stderr info to file path
##          : $stdinfile_path           => Stdinfile path
##          : $stdoutfile_path          => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bedfile;
    my $filehandle;
    my $infile_path;
    my $intermediate_results_dir;
    my $model_type;
    my $num_shards;
    my $outfile_path;
    my $outfile_path_vcf;
    my $referencefile_path;
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
        infile_path => {
            allow       => qr/ bam \z /xms,
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intermediate_results_dir => {
            store       => \$intermediate_results_dir,
            strict_type => 1
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
            store       => \$num_shards,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        outfile_path_vcf => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path_vcf,
            strict_type => 1,
        },
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
      ( get_executable_base_command( { base_command => q{run_deepvariant}, } ), );

    push @commands, q{--reads} . $EQUALS . $infile_path;
    push @commands, q{--ref} . $EQUALS . $referencefile_path;
    push @commands, q{--num_shards} . $EQUALS . $num_shards;
    push @commands, q{--output_gvcf} . $EQUALS . $outfile_path;
    push @commands, q{--output_vcf} . $EQUALS . $outfile_path_vcf;
    push @commands, q{--model_type} . $EQUALS . $model_type;

    if ($bedfile) {
        push @commands, q{--regions} . $EQUALS . $bedfile;
    }
    if ($intermediate_results_dir) {
        push @commands, q{--intermediate_results_dir} . $SPACE . $intermediate_results_dir;
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
