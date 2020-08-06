package MIP::Program::Slurm;

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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ slurm_sacct };
}

sub slurm_sacct {

## Function : Perl wrapper for writing SLURM sacct recipe to already open $filehandle or return commands array. Based on SLURM sacct 2.6.0.
## Returns  : @commands
## Arguments: $fields_format_ref      => List of format fields
##          : $filehandle             => Filehandle to write to
##          : $job_ids_ref            => Slurm job id
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fields_format_ref;
    my $filehandle;
    my $job_ids_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        fields_format_ref => {
            default     => [],
            store       => \$fields_format_ref,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        job_ids_ref => {
            default     => [],
            store       => \$job_ids_ref,
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

    my @commands = q{sacct};

    if ( @{$fields_format_ref} ) {
        push @commands, q{--format=} . join $COMMA, @{$fields_format_ref};
    }

    if ( @{$job_ids_ref} ) {
        push @commands, q{--jobs} . $SPACE . join $COMMA, @{$job_ids_ref};
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

1;
