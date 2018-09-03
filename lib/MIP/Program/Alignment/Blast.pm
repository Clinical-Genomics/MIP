package MIP::Program::Alignment::Blast;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ blast_makeblastdb };
}

## Constants
Readonly my $SPACE => q{ };

sub blast_makeblastdb {

## Function : Perl wrapper for writing blast makeblastdb recipe to $FILEHANDLE. Based on bwa  2.7.1.
## Returns  : @commands
## Arguments: $db_type                => Database type
##          : $FILEHANDLE             => Filehandle to write to
##          : $cdna_seq_file_path     => CDNA sequence file path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $cdna_seq_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $db_type;

    my $tmpl = {
        cdna_seq_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$cdna_seq_file_path,
        },
        db_type => {
            allow       => [qw{ String nucl prot }],
            default     => q{nucl},
            store       => \$db_type,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
