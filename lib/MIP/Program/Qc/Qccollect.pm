package MIP::Program::Qc::Qccollect;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
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
    our @EXPORT_OK = qw{ qccollect };
}

## Constants
Readonly my $SPACE => q{ };

sub qccollect {

## Function : Perl wrapper for qcCollect. Collects information from each analysis run.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $regexp_file_path       => Regular expression file
##          : $log_file_path          => Log file path
##          : $skip_evaluation        => Skip evaluation step

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $infile_path;
    my $outfile_path;
    my $regexp_file_path;
    my $log_file_path;

    ## Default(s)
    my $append_stderr_info;
    my $skip_evaluation;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        regexp_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$regexp_file_path
        },
        stderrfile_path    => { strict_type => 1, store => \$stderrfile_path },
        stdoutfile_path    => { strict_type => 1, store => \$stdoutfile_path },
        log_file_path      => { strict_type => 1, store => \$log_file_path },
        append_stderr_info => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$append_stderr_info
        },
        skip_evaluation => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$skip_evaluation
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{qccollect};

    ## Options
    if ($log_file_path) {

        push @commands, q{--log_file} . $SPACE . $log_file_path;
    }

    if ($regexp_file_path) {

        push @commands, q{--regexp_file} . $SPACE . $regexp_file_path;
    }

    if ($skip_evaluation) {

        push @commands, q{--skip_evaluation} . $SPACE . $skip_evaluation;
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{--sample_info_file} . $SPACE . $infile_path;
    }

    ## Outfile
    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

1;
