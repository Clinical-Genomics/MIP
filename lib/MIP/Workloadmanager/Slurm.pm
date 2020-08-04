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
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      slurm_build_sbatch_header
      slurm_reformat_sacct_output
      slurm_sacct
      slurm_sbatch
    };
}

sub slurm_build_sbatch_header {

## Function : Perl wrapper for writing SLURM sbatch header recipe to already open $filehandle or return commands array. Based on SLURM 2.6.0.
## Returns  : @commands
## Arguments: $core_number              => Core number to allocate
##          : $email                    => User to receive email notification
##          : $email_types_ref          => When to send email for event {REF}
##          : $filehandle               => Filehandle to write to
##          : $job_name                 => Specify a name for the job allocation
##          : $memory_allocation        => Memory allocation
##          : $process_time             => Time limit
##          : $project_id               => Project id
##          : $separator                => Separator to use when writing
##          : $slurm_quality_of_service => Quality of service for the job
##          : $stderrfile_path          => Stderrfile path
##          : $stdoutfile_path          => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $email;
    my $filehandle;
    my $job_name;
    my $memory_allocation;
    my $project_id;
    my $slurm_quality_of_service;
    my $stderrfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $core_number;
    my $process_time;
    my $email_types_ref;
    my $separator;

    my $tmpl = {
        core_number => {
            allow       => qr{ \A\d+\z }xms,
            default     => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        email => {
            store       => \$email,
            strict_type => 1,
        },
        email_types_ref => {
            allow => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref => [qw{ NONE BEGIN END FAIL REQUEUE ALL }],
                            values_ref         => $arg_href->{email_types_ref},
                        }
                    );
                }
            ],
            default     => [qw{ FAIL }],
            store       => \$email_types_ref,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle },
        job_name   => {
            store       => \$job_name,
            strict_type => 1,
        },
        memory_allocation => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            store       => \$memory_allocation,
            strict_type => 1,
        },
        process_time => {
            allow       => [ qr{ \A\d+:\d+:\d+\z }xms, qr{ \A\d-\d+:\d+:\d+\z }xms ],
            default     => q{1:00:00},
            store       => \$process_time,
            strict_type => 1,
        },
        project_id => {
            defined     => 1,
            required    => 1,
            store       => \$project_id,
            strict_type => 1,
        },
        separator => {
            default     => q{\n},
            store       => \$separator,
            strict_type => 1,
        },
        slurm_quality_of_service => {
            store       => \$slurm_quality_of_service,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::List qw{ check_allowed_array_values };
    use MIP::Unix::Write_to_file qw{ unix_write_to_file };

    my @commands;

    push @commands, q{--account=} . $project_id;

    push @commands, q{--ntasks=} . $core_number;

    push @commands, q{--time=} . $process_time;

    if ($memory_allocation) {
        push @commands, q{--mem=} . $memory_allocation . q{G};
    }

    if ($slurm_quality_of_service) {
        push @commands, q{--qos=} . $slurm_quality_of_service;
    }

    if ($job_name) {
        push @commands, q{--job-name=} . $job_name;
    }

    if ($stderrfile_path) {
        push @commands, q{--error=} . $stderrfile_path;
    }

    if ($stdoutfile_path) {
        push @commands, q{--output=} . $stdoutfile_path;
    }

    if ($email) {

        if ( @{$email_types_ref} ) {
            push @commands, q{--mail-type=} . join $COMMA, @{$email_types_ref};
        }
        push @commands, q{--mail-user=} . $email;
    }

    ## Add sbatch shebang
    @commands = map { q{#SBATCH } . $_ } @commands;

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $separator,
        }
    );
    return @commands;
}

sub slurm_sacct {

## Function : Perl wrapper for writing SLURM sacct recipe to already open $filehandle or return commands array. Based on SLURM sacct 2.6.0.
## Returns  : "@commands"
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

    ## sacct
    my @commands = q{sacct};    #Stores commands depending on input parameters

    ## Options
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

sub slurm_reformat_sacct_output {

## Function : Removes ".batch" lines in sacct output.
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

    #Stream to
    print {$filehandle} $TAB . join( $SPACE, @{$commands_ref} ) . $SPACE;

    # Pipe
    print {$filehandle} $PIPE . $SPACE;

    #Execute perl
    print {$filehandle} q?perl -nae '?;

    # Define reformated sacct headers
    print {$filehandle} q?my @headers=(?
      . join( $COMMA, @{$reformat_sacct_headers_ref} ) . q?); ?;

    #Write header line
    print {$filehandle} q?if($. == 1) { print q{#} . join(qq{\t}, @headers), qq{\n} } ?;

  # Write individual job line - skip line containing (.batch or .bat+) in the first column
    print {$filehandle}
q?if ($. >= 3 && $F[0] !~ /( .batch | .bat+ )\b/xms) { print join(qq{\t}, @F), qq{\n} }' ?;

    #Write to log_file.status
    print {$filehandle} q{> } . $log_file_path . q{.status} . $NEWLINE x 2;
    return;
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

1;
