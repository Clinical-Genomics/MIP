package MIP::Program::Variantcalling::Star_fusion;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ star_fusion };
}

## Constants
Readonly my $SPACE => q{ };

sub star_fusion {

## Function : Splice site analysis using STAR-fusion.
## Returns  :
## Arguments: $cpu                    => Number of threads for running STAR
##          : $examine_coding_effect  => Append coding effect to fusion
##          : $fastq_r1_path          => The path of the R1 fastq
##          : $fastq_r2_path          => The path of the R2 fastq
##          : $FILEHANDLE             => Filehandle to write to
##          : $genome_lib_dir_path    => Path to the directory containing the genome library
##          : $output_directory_path  => output directory path
##          : $samples_file_path      => Sample file path
##          : $sjdb_path              => Splice junction database file path (the junctions tab file)
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cpu;
    my $fastq_r1_path;
    my $fastq_r2_path;
    my $FILEHANDLE;
    my $genome_lib_dir_path;
    my $output_directory_path;
    my $samples_file_path;
    my $sjdb_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $examine_coding_effect;

    my $tmpl = {
        cpu => {
            store       => \$cpu,
            strict_type => 1,
        },
        examine_coding_effect => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$examine_coding_effect,
            strict_type => 1,
        },
        fastq_r1_path => {
            store       => \$fastq_r1_path,
            strict_type => 1,
        },
        fastq_r2_path => {
            store       => \$fastq_r2_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        genome_lib_dir_path => {
            required    => 1,
            store       => \$genome_lib_dir_path,
            strict_type => 1,
        },
        output_directory_path => {
            required    => 1,
            store       => \$output_directory_path,
            strict_type => 1,
        },
        samples_file_path => {
            store       => \$samples_file_path,
            strict_type => 1,
        },
        sjdb_path => {
            store       => \$sjdb_path,
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
    my @commands = q{STAR-Fusion};

    push @commands, q{--genome_lib_dir} . $SPACE . $genome_lib_dir_path;

    if ($samples_file_path) {

        push @commands, q{--samples_file} . $SPACE . $samples_file_path;
    }
    elsif ( $sjdb_path and not $fastq_r1_path ) {

        push @commands, q{-J} . $SPACE . $sjdb_path;
    }
    elsif ( $fastq_r1_path and $fastq_r2_path ) {

        push @commands, q{--right_fq} . $SPACE . $fastq_r1_path;
        push @commands, q{--left_fq} . $SPACE . $fastq_r2_path;
    }
    else {
        croak(
q{Error: You must either specify the fastq file paths or a splice junction database (SJDB) file.}
        );
    }

    if ($cpu) {
        push @commands, q{--CPU} . $SPACE . $cpu;
    }

    if ($examine_coding_effect) {
        push @commands, q{--examine_coding_effect};
    }

    push @commands, q{--output_dir} . $SPACE . $output_directory_path;

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
