package MIP::Program::Vcf_compare;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ vcf_compare };
}

sub vcf_compare {

## Function : Perl wrapper for generic commands module.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $infile_path_truth      => Truth infile vcf
##          : $infile_path_vcf        => vcf-file to compare with truth set

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $infile_path_truth;
    my $infile_path_vcf;

    ## Default(s)

    my $tmpl = {
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
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        infile_path_truth => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path_truth,
            strict_type => 1,
        },
        infile_path_vcf => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path_vcf,
            strict_type => 1,
        }


    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ vcf_compare };

    push @commands, $infile_path_truth;
    push @commands, $infile_path_vcf;

    ############################################
    ## ADD COMMAND SPECIFIC FLAGS AND OPTIONS ##
    ############################################

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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
