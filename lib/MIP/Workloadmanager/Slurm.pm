package MIP::Workloadmanager::Slurm;

#### Copyright 2017 Henrik Stranneheim

use strict;
use warnings;
use warnings qw(FATAL utf8);
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

use FindBin qw($Bin);                 #Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw(unix_standard_streams);
use MIP::Unix::Write_to_file qw(unix_write_to_file);

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(slurm_sacct slurm_sbatch slurm_build_sbatch_header);
}

sub slurm_sacct {

##slurm_sacct

##Function : Perl wrapper for writing SLURM sacct recipe to already open $FILEHANDLE or return commands array. Based on SLURM sacct 2.6.0.
##Returns  : "@commands"
##Arguments: $fields_format_ref, $job_ids_ref, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $stderrfile_path_append,
##         : $job_ids_ref            => Slurm job id
##         : $fields_format_ref      => List of format fields
##         : $stderrfile_path        => Stderrfile path
##         : $stdoutfile_path        => Stdoutfile path
##         : $FILEHANDLE             => Filehandle to write to
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fields_format_ref;
    my $job_ids_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        fields_format_ref =>
          { default => [], strict_type => 1, store => \$fields_format_ref },
        job_ids_ref =>
          { default => [], strict_type => 1, store => \$job_ids_ref },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $SPACE = q{ };
    my $COMMA = q{,};

    ## sacct
    my @commands = qw(sacct);    #Stores commands depending on input parameters

    ## Options
    if ( @{$fields_format_ref} ) {

        push @commands, '--format=' . join $COMMA, @{$fields_format_ref};
    }
    if ( @{$job_ids_ref} ) {

        push @commands, '--jobs ' . join $COMMA, @{$job_ids_ref};
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

sub slurm_sbatch {

##slurm_sbatch

##Function : Perl wrapper for writing SLURM sbatch recipe to already open $FILEHANDLE or return commands array. Based on SLURM sbatch 2.6.0.
##Returns  : "@commands"
##Arguments: infile_path, $dependency_type, $job_ids_string, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $stderrfile_path_append,
##         : $infile_path            => Infile_path
##         : $dependency_type        => Type of slurm job dependency
##         : $job_ids_string         => Slurm job ids string (:job_id:job_id...n)
##         : $stderrfile_path        => Stderrfile path
##         : $stdoutfile_path        => Stdoutfile path
##         : $FILEHANDLE             => Filehandle to write to
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Default(s)
    my $stderrfile_path_append;

    ## Flatten argument(s)
    my $infile_path;
    my $dependency_type;
    my $job_ids_string;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        dependency_type => { strict_type => 1, store => \$dependency_type },
        job_ids_string  => { strict_type => 1, store => \$job_ids_string },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $SPACE = q{ };

    ## sbatch
    my @commands = qw(sbatch);    #Stores commands depending on input parameters

    ## Options
    if ( ($dependency_type) && ($job_ids_string) ) {

        push @commands, '--dependency=' . $dependency_type . $job_ids_string;
    }

    # Infile
    push @commands, $infile_path;

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

sub slurm_build_sbatch_header {

##slurm_build_sbatch_header

##Function : Perl wrapper for writing SLURM sbatch header recipe to already open $FILEHANDLE or return commands array. Based on SLURM 2.6.0.
##Returns  : "@commands"
##Arguments: $project_id, $slurm_quality_of_service, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $core_number, $process_time, $job_name, $email, $email_types_ref, $separator
##         : $project_id               => Project id
##         : $slurm_quality_of_service => Quality of service for the job
##         : $job_name                 => Specify a name for the job allocation
##         : $stderrfile_path          => Stderrfile path
##         : $stdoutfile_path          => Stdoutfile path
##         : $email                    => User to receive email notification
##         : $FILEHANDLE               => Filehandle to write to
##         : $core_number              => Core number to allocate
##         : $process_time             => Time limit
##         : $email_types_ref          => When to send email for event
##         : $separator                => Separator to use when writing

    my ($arg_href) = @_;

    ## Default(s)
    my $core_number;
    my $process_time;
    my $email_types_ref;
    my $separator;

    ## Flatten argument(s)
    my $project_id;
    my $slurm_quality_of_service;
    my $job_name;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $email;
    my $FILEHANDLE;

    ## Constants
    my $NEWLINE = q{\n};

    use MIP::Check::Parameter qw(check_allowed_array_values);

    my $tmpl = {
        project_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$project_id
        },
        slurm_quality_of_service =>
          { strict_type => 1, store => \$slurm_quality_of_service },
        job_name        => { strict_type => 1, store => \$job_name },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        email           => { strict_type => 1, store => \$email },
        email_types_ref => {
            default => ['FAIL'],
            allow   => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref =>
                              [qw(NONE BEGIN END FAIL REQUEUE ALL)],
                            values_ref => $arg_href->{email_types_ref},
                        }
                    );
                }
            ],
            strict_type => 1,
            store       => \$email_types_ref
        },
        FILEHANDLE  => { store => \$FILEHANDLE },
        core_number => {
            default     => 1,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$core_number
        },
        process_time => {
            default     => '1:00:00',
            allow       => qr/^\d+:\d+:\d+$/,
            strict_type => 1,
            store       => \$process_time
        },
        separator => {
            default     => $NEWLINE,
            strict_type => 1,
            store       => \$separator
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    use MIP::Unix::Write_to_file qw(unix_write_to_file);

    ## Build sbatch header
    my @commands;    #Stores commands depending on input parameters

    ## Options
    # Account allocation
    push @commands, '--account=' . $project_id;

    # Core allocation
    push @commands, '--ntasks=' . $core_number;

    # Time allocation
    push @commands, '--time=' . $process_time;

    # Quality of service
    if ($slurm_quality_of_service) {

        push @commands, '--qos=' . $slurm_quality_of_service;
    }

    # Job name
    if ($job_name) {

        push @commands, '--job-name=' . $job_name;
    }

    # Stderror
    if ($stderrfile_path) {

        push @commands, '--error=' . $stderrfile_path;
    }

    # Stdout
    if ($stdoutfile_path) {

        push @commands, '--output=' . $stdoutfile_path;
    }

    # Email and email_types
    if ($email) {

        if ( @{$email_types_ref} ) {

            push @commands, '--mail-type=' . join q{,}, @{$email_types_ref};
        }
        push @commands, '--mail-user=' . $email;
    }

    ## Add sbatch shebang
    @commands = map { '#SBATCH ' . $_ } @commands;

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $separator,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

1;
