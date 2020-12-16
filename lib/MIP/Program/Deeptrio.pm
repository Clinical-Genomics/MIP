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
use MIP::Constants qw{ $SPACE $EQUALS };
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
##          : $iofile_parameters_href  => Path to input and output files, and sample ids
##          : $model_type             => Type of model to use for variant calling. Allowed values WES, WGS, or PACBIO
##          : $num_shards             => Number of files the input is split into for the make examples step
##          : $referencefile_path     => Path to genome reference
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bedfile;
    my $filehandle;
    my $iofile_parameters_href;
    my $model_type;
    my $num_shards;
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
        iofile_parameters_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$iofile_parameters_href,
            strict_type => 1,
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
      ( get_executable_base_command( { base_command => q{run_deeptrio}, } ), );

    push @commands, q{--model_type} . $SPACE . $model_type;
    push @commands, q{--ref} . $SPACE . $referencefile_path;
    push @commands, q{--num_shards} . $SPACE . $num_shards;
    push @commands, q{--reads_child} . $SPACE . $iofile_parameters_href->{reads_child};
    push @commands, q{--sample_name_child} . $SPACE . $iofile_parameters_href->{sample_name_child};
    push @commands, q{--output_gvcf_child} . $SPACE . $iofile_parameters_href->{output_gvcf_child};
    push @commands, q{--output_vcf_child} . $SPACE . $iofile_parameters_href->{output_vcf_child};
    push @commands, q{--reads_parent1} . $SPACE . $iofile_parameters_href->{reads_parent1};
    push @commands, q{--sample_name_parent1} . $SPACE . $iofile_parameters_href->{sample_name_parent1};
    push @commands, q{--output_gvcf_parent1} . $SPACE . $iofile_parameters_href->{output_gvcf_parent1};
    push @commands, q{--output_vcf_parent1} . $SPACE . $iofile_parameters_href->{output_vcf_parent1};
    if ($iofile_parameters_href->{sample_name_parent2}){
        push @commands, q{--reads_parent2} . $SPACE . $iofile_parameters_href->{reads_parent2};
        push @commands, q{--sample_name_parent2} . $SPACE . $iofile_parameters_href->{sample_name_parent2};
        push @commands, q{--output_gvcf_parent2} . $SPACE . $iofile_parameters_href->{output_gvcf_parent2};
        push @commands, q{--output_vcf_parent2} . $SPACE . $iofile_parameters_href->{output_vcf_parent2};
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
