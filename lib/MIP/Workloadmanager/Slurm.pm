package MIP::Workloadmanager::Slurm;

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
use MIP::Constants qw{ $COMMA $EQUALS $NEWLINE $PIPE $SPACE $TAB };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.11;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      slurm_sbatch
      slurm_reformat_sacct_output
      slurm_track_progress
    };
}

sub slurm_sbatch {

## Function : Perl wrapper for writing SLURM sbatch recipe to already open $filehandle or return commands array. Based on SLURM sbatch 2.6.0.
## Returns  : @commands
## Arguments: $base_command           => Sbatch or slurm-mock (for tests)
##          : $dependency_type        => Type of slurm job dependency
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Infile_path
##          : $job_ids_string         => Slurm job ids string (:job_id:job_id...n)
##          : $reservation_name       => Allocate resources from named reservation
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dependency_type;
    my $filehandle;
    my $infile_path;
    my $job_ids_string;
    my $reservation_name;
    my $stderrfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $base_command;
    my $stderrfile_path_append;

    my $tmpl = {
        base_command => {
            default     => q{sbatch},
            store       => \$base_command,
            strict_type => 1,
        },
        dependency_type => {
            store       => \$dependency_type,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        job_ids_string => {
            store       => \$job_ids_string,
            strict_type => 1,
        },
        reservation_name => {
            store       => \$reservation_name,
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

    my @commands = $base_command;

    if ( $dependency_type and $job_ids_string ) {

        push @commands, q{--dependency=} . $dependency_type . $job_ids_string;
    }

    if ($reservation_name) {

        push @commands, q{--reservation} . $EQUALS . $reservation_name;
    }
    push @commands, $infile_path;

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

sub slurm_reformat_sacct_output {

## Function : Removes ".batch" lines in sacct output
## Returns  :
## Arguments: $commands_ref               => Commands to stream to perl oneliner
##          : $filehandle                 => Sbatch filehandle to write to
##          : $log_file_path              => The log file {REF}
##          : $reformat_sacct_headers_ref => Reformated sacct headers

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $filehandle;
    my $log_file_path;
    my $reformat_sacct_headers_ref;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        log_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$log_file_path,
            strict_type => 1,
        },
        reformat_sacct_headers_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$reformat_sacct_headers_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };

    # Stream to
    print {$filehandle} $TAB . join( $SPACE, @{$commands_ref} ) . $SPACE;

    # Pipe
    print {$filehandle} $PIPE . $SPACE;

    perl_nae_oneliners(
        {
            filehandle         => $filehandle,
            oneliner_name      => q{reformat_sacct_headers},
            oneliner_parameter => join( $COMMA, @{$reformat_sacct_headers_ref} ),
            stdoutfile_path    => $log_file_path . q{.status} . $NEWLINE x 2,
        }
    );
    return;
}

sub slurm_track_progress {

## Function : Output Slurm info on each job via sacct command and write to log file path
## Returns  :
## Arguments: $filehandle              => Sbatch filehandle to write to
##          : $job_ids_ref             => Job ids
##          : $log_file_path           => The log file path to write status to
##          : $sacct_format_fields_ref => Format and fields of sacct output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $job_ids_ref;
    my $log_file_path;

    ## Default(s)
    my $sacct_format_fields_ref;

    my $tmpl = {
        filehandle    => { store   => \$filehandle, },
        job_ids_ref   => { default => [], store => \$job_ids_ref, strict_type => 1, },
        log_file_path => { store   => \$log_file_path, strict_type => 1, },
        sacct_format_fields_ref => {
            default => [
                qw{
                  jobid jobname%50 account partition alloccpus TotalCPU elapsed start end state exitcode }
            ],
            store       => \$sacct_format_fields_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Slurm qw{ slurm_sacct };
    use MIP::Workloadmanager::Slurm qw{ slurm_reformat_sacct_output };

    return if ( not @{$job_ids_ref} );

    my @reformat_sacct_headers = @{$sacct_format_fields_ref};

  HEADER_ELEMENT:
    foreach my $element (@reformat_sacct_headers) {

        ## Remove "%digits" from headers
        $element =~ s/%\d+//gsxm;
    }

    my @commands = slurm_sacct(
        {
            fields_format_ref => \@{$sacct_format_fields_ref},
            job_ids_ref       => \@{$job_ids_ref},
        }
    );

    slurm_reformat_sacct_output(
        {
            commands_ref               => \@commands,
            filehandle                 => $filehandle,
            log_file_path              => $log_file_path,
            reformat_sacct_headers_ref => \@reformat_sacct_headers,
        }
    );
    return;
}

1;
