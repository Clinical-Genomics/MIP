package MIP::Program::Blast;

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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ blast_blastn blast_makeblastdb };
}

sub blast_blastn {

## Function : Perl wrapper for writing blast blastn recipe to $filehandle. Based on bwa  2.7.1.
## Returns  : @commands
## Arguments: $database_name          => Database type
##          : $evalue                 => Expectation value (E) threshold for saving hits
##          : $filehandle             => Filehandle to write to
##          : $lcase_masking          => Use lower case filtering in query and subject sequence(s)
##          : $max_target_seqs        => Max number of target sequences
##          : $output_format          => Output format
##          : $query_file_path        => Query file path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $thread_number          => Number of threads (CPUs) to use in the BLAST search
##          : $word_size              => Word size for wordfinder algorithm (length of best perfect match

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $query_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $database_name;
    my $evalue;
    my $lcase_masking;
    my $max_target_seqs;
    my $output_format;
    my $thread_number;
    my $word_size;

    ## Constants
    Readonly my $EXPECT_VALUE      => 1e-3;
    Readonly my $MAX_TARGET_SEQS   => 1000;
    Readonly my $MAX_THREAD_NUMBER => 36;
    Readonly my $TABULAR           => 6;
    Readonly my $THREAD_NUMBER     => 1;
    Readonly my $WORD_SIZE         => 11;

    my $tmpl = {
        evalue => {
            default     => $EXPECT_VALUE,
            store       => \$evalue,
            strict_type => 1,
        },
        database_name => {
            defined     => 1,
            required    => 1,
            store       => \$database_name,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        lcase_masking => {
            allow       => [ 0, 1, undef ],
            default     => 0,
            store       => \$lcase_masking,
            strict_type => 1,
        },
        max_target_seqs => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $MAX_TARGET_SEQS,
            store       => \$max_target_seqs,
            strict_type => 1,
        },
        query_file_path => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$query_file_path,
        },
        output_format => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $TABULAR,
            store       => \$output_format,
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
            allow       => [ 1 .. $MAX_THREAD_NUMBER ],
            default     => $THREAD_NUMBER,
            store       => \$thread_number,
            strict_type => 1,
        },
        word_size => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $WORD_SIZE,
            store       => \$word_size,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ blastn };

    push @commands, q{-query} . $SPACE . $query_file_path;

    push @commands, q{-db} . $SPACE . $database_name;

    push @commands, q{-max_target_seqs} . $SPACE . $max_target_seqs;

    push @commands, q{-outfmt} . $SPACE . $output_format;

    push @commands, q{-evalue} . $SPACE . $evalue;

    push @commands, q{-num_threads} . $SPACE . $thread_number;

    push @commands, q{-word_size} . $SPACE . $word_size;

    if ($lcase_masking) {

        push @commands, q{-lcase_masking};
    }

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
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub blast_makeblastdb {

## Function : Perl wrapper for writing blast makeblastdb recipe to $filehandle. Based on bwa  2.7.1.
## Returns  : @commands
## Arguments: $db_type                => Database type
##          : $filehandle             => Filehandle to write to
##          : $cdna_seq_file_path     => CDNA sequence file path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $cdna_seq_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $db_type;

    my $tmpl = {
        cdna_seq_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$cdna_seq_file_path,
            strict_type => 1,
        },
        db_type => {
            allow       => [qw{ String nucl prot }],
            default     => q{nucl},
            store       => \$db_type,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
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
    my @commands = qw{ makeblastdb };

    push @commands, q{-in} . $SPACE . $cdna_seq_file_path;

    push @commands, q{-dbtype} . $SPACE . $db_type;

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
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
