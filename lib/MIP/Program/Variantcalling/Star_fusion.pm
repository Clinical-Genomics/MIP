package MIP::Program::Variantcalling::Star_fusion;

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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ star_fusion };
}

## Constants
Readonly my $SPACE => q{ };

sub star_fusion {

## Function :
## Returns  :
## Arguments: $fastq_r1_path          => The path of the R1 fastq
##          : $fastq_r2_path          => The path of the R2 fastq
##          : $FILEHANDLE             => Filehandle to write to
##          : $genome_lib_dir_path    => Path to the directory containing the genome library
##          : $output_directory_path  => output directory path
##          : $sjdb_path              => Splice junction database file path (the junctions tab file)
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fastq_r1_path;
    my $fastq_r2_path;
    my $FILEHANDLE;
    my $genome_lib_dir_path;
    my $sjdb_path;
    my $output_directory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        fastq_r1_path  => { strict_type => 1, store => \$fastq_r1_path, },
        fastq_r2_path => { strict_type => 1, store => \$fastq_r2_path, },
        FILEHANDLE     => {
            store => \$FILEHANDLE,
        },
        genome_lib_dir_path =>
          { required => 1, strict_type => 1, store => \$genome_lib_dir_path, },
        output_directory_path => {
            required    => 1,
            strict_type => 1,
            store       => \$output_directory_path,
        },
        sjdb_path       => { strict_type => 1, store => \$sjdb_path, },
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{STAR-Fusion};

    push @commands, q{--genome_lib_dir} . $SPACE . $genome_lib_dir_path;
    if ( $sjdb_path and not $fastq_r1_path ) {
        push @commands, q{-J} . $SPACE . $sjdb_path;
    }
    elsif ( $fastq_r1_path and $fastq_r2_path ) {
        push @commands,
            q{--right_fq}
          . $SPACE
          . $fastq_r1_path
          . $SPACE
          . q{--left_fq}
          . $SPACE
          . $fastq_r2_path;
    }
    else {
        croak
q{Error: You must either specify the fastq file paths or a splice junction database (SJDB) file.};
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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
