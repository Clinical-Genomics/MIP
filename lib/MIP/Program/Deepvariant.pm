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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ deepvariant };
}

sub deepvariant {

## Function : Perl wrapper for generic commands module
## Returns  : @commands
## Arguments: $bamfile                => Aligned, sorted, indexed bam file containing the reads we want to call
##          : $bedfile                => Bed file containing the list of  regions we want to process
##          : $filehandle             => Filehandle to write to
##          : $model_type             => Type of model to use for variant calling. Allowed values WES, WGS, or PACBIO
##          : $num_shards             => Number of files the input is split into for the make examples step
##          : $outfile_path           => Path to the output gvcf file
##          : $referencefile_path     => Path to genome reference
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bamfile;
    my $bedfile;
    my $filehandle;
    my $model_type;
    my $num_shards;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    my $tmpl = {
        bamfile => {
            allow       => qr/ bam \z /xms,
            defined     => 1,
            required    => 1,
            store       => \$bamfile,
            strict_type => 1,
        },
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
            store       => \$num_shards,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
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

    ## Stores commands depending on input parameters
    my @commands = qw{ run_deepvariant };

    push @commands, q{--reads} . $EQUALS . $bamfile;
    push @commands, q{--ref} . $EQUALS . $referencefile_path;
    push @commands, q{--output_gvcf} . $EQUALS . $outfile_path;
    push @commands, q{--model_type} . $EQUALS . $model_type;

    if ($bedfile) {
        push @commands, q{--regions} . $EQUALS . $bedfile;
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
