package File::Format::Shell;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use v5.10;    #Require at least perl 5.10
use utf8;     #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.01;

    # Inherit from Exporter to export functions and variables
    our @ISA = qw(Exporter);

    # Functions and variables which are exported by default
    our @EXPORT = qw();

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = (
        'create_bash_file',             'build_shebang',
        'create_housekeeping_function', 'create_error_trap_function',
        'enable_trap',                  'clear_trap',
        'track_progress',
    );
}

use Cwd;
use FindBin qw($Bin);    #Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catfile catdir devnull);
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

##MIPs lib/
use lib catdir( $Bin, 'lib' );        #Add MIPs internal lib

sub create_bash_file {

##create_bash_file

##Function : Create bash file with header
##Returns  : "$FILEHANDLE"
##Arguments: $file_name, $directory_remove, $log, $trap_signals_ref, $trap_function, $set_login_shell, $set_errexit, $set_nounset, $set_pipefail
##         : $file_name        => File name
##         : $directory_remove => Directory to remove when caught by trap function
##         : $log              => Log object to write to
##         : $trap_signals_ref => Array with signals to clear trap for {REF}
##         : $trap_function    => Trap function argument
##         : $set_login_shell  => Invoked as a login shell. Reinitilize bashrc and bash_profile
##         : $set_errexit      => Halt script if command has non-zero exit code (-e)
##         : $set_nounset      => Halt script if variable is uninitialised (-u)
##         : $set_pipefail     => Detect errors within pipes (-o pipefail)

    my ($arg_href) = @_;

    ## Default(s)
    my $trap_signals_ref;
    my $trap_function;
    my $set_login_shell;
    my $set_errexit;
    my $set_nounset;
    my $set_pipefail;

    ## Flatten argument(s)
    my $file_name;
    my $directory_remove;
    my $log;

    my $tmpl = {
        file_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_name
        },
        log              => { store => \$log },
        directory_remove => {
            allow       => qr/^\.\S+$/,
            strict_type => 1,
            store       => \$directory_remove
        },
        trap_signals_ref => {
            default     => ['ERR'],
            strict_type => 1,
            store       => \$trap_signals_ref
        },
        trap_function => {
            default     => 'error',
            strict_type => 1,
            store       => \$trap_function
        },
        set_login_shell => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_login_shell
        },
        set_errexit => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_errexit
        },
        set_nounset => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_nounset
        },
        set_pipefail => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_pipefail
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    ## Open batch file with supplied log
    if ( ( defined $log ) && ($log) ) {

        open( $FILEHANDLE, '>', catfile($file_name) )
          or $log->logdie(
            "Cannot write to '" . catfile($file_name) . "' :" . $! . "\n" );
    }
    else {

        open( $FILEHANDLE, '>', catfile($file_name) )
          or
          die( "Cannot write to '" . catfile($file_name) . "' :" . $! . "\n" );
    }

    # Build bash shebang line
    build_shebang(
        {
            FILEHANDLE      => $FILEHANDLE,
            set_login_shell => $set_login_shell,
            set_errexit     => $set_errexit,
            set_nounset     => $set_nounset,
            set_pipefail    => $set_pipefail,
        }
    );

    ## Create housekeeping function which removes entire directory when finished
    create_housekeeping_function(
        {
            directory_remove => $directory_remove,
            FILEHANDLE       => $FILEHANDLE,
        }
    );

    ## Create debug trap
    enable_trap(
        {
            FILEHANDLE       => $FILEHANDLE,
            trap_signals_ref => ["DEBUG"],
            trap_function    => q?previous_command="$BASH_COMMAND"?,
        }
    );

    ## Create error handling function and trap
    create_error_trap_function( { FILEHANDLE => $FILEHANDLE, } );

    if ( ( defined($log) ) && ($log) ) {

        $log->info( "Created bash file: '" . catfile($file_name), "'\n" );
    }
    else {

        print STDERR "Created bash file: '" . catfile($file_name), "'", "\n";
    }
    return $FILEHANDLE;
}

sub build_shebang {

##build_shebang

##Function : Build bash shebang line
##Returns  : ""
##Arguments: $FILEHANDLE, $bash_bin_path, $set_login_shell, $set_errexit, $set_nounset, $set_pipefail
##         : $FILEHANDLE      => Filehandle to write to
##         : $bash_bin_path   => Location of bash bin
##         : $set_login_shell => Invoked as a login shell (-l). Reinitilize bashrc and bash_profile
##         : $set_errexit     => Halt script if command has non-zero exit code (-e)
##         : $set_nounset     => Halt script if variable is uninitialised (-u)
##         : $set_pipefail    => Detect errors within pipes (-o pipefail)

    my ($arg_href) = @_;

    ## Default(s)
    my $bash_bin_path;
    my $set_login_shell;
    my $set_errexit;
    my $set_nounset;
    my $set_pipefail;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
        FILEHANDLE    => { required => 1, store => \$FILEHANDLE },
        bash_bin_path => {
            default =>
              catfile( dirname( dirname( devnull() ) ), qw(usr bin env bash) ),
            allow       => qr/^\S+$/,
            strict_type => 1,
            store       => \$bash_bin_path
        },
        set_login_shell => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_login_shell
        },
        set_errexit => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_errexit
        },
        set_nounset => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_nounset
        },
        set_pipefail => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_pipefail
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    # Build shebang
    print $FILEHANDLE '#!';
    print $FILEHANDLE $bash_bin_path;
    print $FILEHANDLE "\n";

    # Set flags
    print $FILEHANDLE "set -l\n"          if ($set_login_shell);
    print $FILEHANDLE "set -e\n"          if ($set_errexit);
    print $FILEHANDLE "set -u\n"          if ($set_nounset);
    print $FILEHANDLE "set -o pipefail\n" if ($set_pipefail);
    print $FILEHANDLE "\n";
    return;
}

sub create_housekeeping_function {

##create_housekeeping_function

##Function : Create housekeeping function which removes entire directory when finished
##Returns  : ""
##Arguments: $job_ids_ref, $sacct_format_fields_ref, $log_file_ref, $FILEHANDLE, $directory_remove, $trap_signals_ref, $trap_function
##         : $job_ids_ref             => Job ids
##         : $sacct_format_fields_ref => Format and fields of sacct output
##         : $log_file_ref            => Log file to write job_id progress to {REF}
##         : $FILEHANDLE              => Filehandle to write to
##         : $directory_remove        => Directory to remove when caught by trap function
##         : $trap_signals_ref        => Array with signals to enable trap for {REF}
##         : $trap_function           => The trap function argument

    my ($arg_href) = @_;

    ## Default(s)
    my $trap_signals_ref;
    my $trap_function;

    ## Flatten argument(s)
    my $job_ids_ref;
    my $sacct_format_fields_ref;
    my $log_file_ref;
    my $FILEHANDLE;
    my $directory_remove;

    my $tmpl = {
        job_ids_ref =>
          { default => [], strict_type => 1, store => \$job_ids_ref },
        sacct_format_fields_ref => {
            default     => [],
            strict_type => 1,
            store       => \$sacct_format_fields_ref
        },
        log_file_ref =>
          { default => \$$, strict_type => 1, store => \$log_file_ref },
        FILEHANDLE       => { required    => 1, store => \$FILEHANDLE },
        directory_remove => { strict_type => 1, store => \$directory_remove },
        trap_signals_ref => {
            default     => [ "EXIT", "TERM", "INT" ],
            strict_type => 1,
            store       => \$trap_signals_ref
        },
        trap_function => {
            default     => "finish",
            strict_type => 1,
            store       => \$trap_function
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    use Program::Gnu::Coreutils qw(rm);

    ## Create housekeeping function and trap
    print $FILEHANDLE q?finish() {?, "\n\n";

    if ( ( defined($directory_remove) ) && ($directory_remove) ) {

        print $FILEHANDLE "\t" . q?local directory="$1"?,         "\n";
        print $FILEHANDLE "\t" . q?## Perform exit housekeeping?, "\n";
        print $FILEHANDLE "\t";
        rm(
            {
                infile_path => q?"$directory"?,
                force       => 1,
                recursive   => 1,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        print $FILEHANDLE "\n\n";
    }
    if (   ( defined($job_ids_ref) )
        && (@$job_ids_ref)
        && ( defined($$log_file_ref) )
        && ($$log_file_ref) )
    {

        ## Output SLURM info on each job via sacct command and write to log file(.status)
        track_progress(
            {
                job_ids_ref             => \@{$job_ids_ref},
                sacct_format_fields_ref => \@{$sacct_format_fields_ref},
                FILEHANDLE              => $FILEHANDLE,
                log_file_ref            => $log_file_ref,
            }
        );
    }

    print $FILEHANDLE q?}?, "\n";

    ## Enable trap function with trap signal(s)
    enable_trap(
        {
            FILEHANDLE       => $FILEHANDLE,
            trap_signals_ref => \@{$trap_signals_ref},
            trap_function    => $trap_function,
        }
    );
}

sub create_error_trap_function {

##create_error_trap_function

##Function : Create error handling function and trap
##Returns  : ""
##Arguments: $job_ids_ref, sacct_format_fields_ref, $log_file_ref, $FILEHANDLE, $trap_function_call, $trap_signals_ref, $trap_function_name
##         : $job_ids_ref             => Job ids
##         : $sacct_format_fields_ref => Format and fields of sacct output
##         : $log_file_ref            => Log file to write job_id progress to {REF}
##         : $FILEHANDLE              => Filehandle to write to
##         : $trap_function_call      => Trap function call
##         : $trap_signals_ref        => Array with signals to enable trap for {REF}
##         : $trap_function_name      => The trap function argument

    my ($arg_href) = @_;

    ## Default(s)
    my $trap_signals_ref;
    my $trap_function_name;

    ## Flatten argument(s)
    my $job_ids_ref;
    my $sacct_format_fields_ref;
    my $log_file_ref;
    my $FILEHANDLE;
    my $trap_function_call;

    my $tmpl = {
        job_ids_ref =>
          { default => [], strict_type => 1, store => \$job_ids_ref },
        sacct_format_fields_ref => {
            default     => [],
            strict_type => 1,
            store       => \$sacct_format_fields_ref
        },
        log_file_ref =>
          { default => \$$, strict_type => 1, store => \$log_file_ref },
        FILEHANDLE         => { required => 1, store => \$FILEHANDLE },
        trap_function_call => {
            default     => q{$(error "$previous_command" "$?")},
            strict_type => 1,
            store       => \$trap_function_call
        },
        trap_signals_ref => {
            default     => ["ERR"],
            strict_type => 1,
            store       => \$trap_signals_ref
        },
        trap_function_name => {
            default     => "error",
            strict_type => 1,
            store       => \$trap_function_name
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    ## Create error handling function and trap
    print $FILEHANDLE $trap_function_name . q?() {?, "\n\n";
    print $FILEHANDLE "\t" . q?local program="$1"?,     "\n";
    print $FILEHANDLE "\t" . q?local return_code="$2"?, "\n\n";

    if (   ( defined($job_ids_ref) )
        && (@$job_ids_ref)
        && ( defined($$log_file_ref) )
        && ($$log_file_ref) )
    {

        ## Output SLURM info on each job via sacct command and write to log file(.status)
        track_progress(
            {
                job_ids_ref             => \@{$job_ids_ref},
                sacct_format_fields_ref => \@{$sacct_format_fields_ref},
                FILEHANDLE              => $FILEHANDLE,
                log_file_ref            => $log_file_ref,
            }
        );
    }

    print $FILEHANDLE "\t" . q?## Display error message and exit?, "\n";
    print $FILEHANDLE "\t"
      . q?echo "${program}: ${return_code}: Unknown Error - ExitCode=$return_code" 1>&2?,
      "\n";
    print $FILEHANDLE "\t" . q?exit 1?, "\n";
    print $FILEHANDLE q?}?, "\n";

    ## Enable trap function with trap signal(s)
    enable_trap(
        {
            FILEHANDLE       => $FILEHANDLE,
            trap_signals_ref => \@{$trap_signals_ref},
            trap_function    => $trap_function_call,
        }
    );
}

sub clear_trap {

##clear_trap

##Function : Clear trap for signal(s), e.g. in exome analysis since the might be no variants in MT or Y contigs. This will cause premature exit from sbatch
##Returns  : ""
##Arguments: $FILEHANDLE, $trap_signals_ref
##         : $FILEHANDLE       => The FILEHANDLE to write to
##         : $trap_signals_ref => Array with signals to clear trap for {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $trap_signals_ref;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
        FILEHANDLE       => { required => 1, store => \$FILEHANDLE },
        trap_signals_ref => {
            default     => ["ERR"],
            strict_type => 1,
            store       => \$trap_signals_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    ## Clear trap for signal ERR
    print $FILEHANDLE "\n## Clear trap for signal(s) "
      . join( " ", @$trap_signals_ref ), "\n";
    print $FILEHANDLE "trap - " . join( " ", @$trap_signals_ref ), "\n";
    print $FILEHANDLE "trap", "\n\n";
}

sub enable_trap {

##enable_trap

##Function : Enable trap function with trap signal(s).
##Returns  : ""
##Arguments: $FILEHANDLE, $trap_signals_ref, $trap_function
##         : $FILEHANDLE       => The FILEHANDLE to write to
##         : $trap_signals_ref => Array with signals to enable trap for {REF}
##         : $trap_function    => The trap function argument

    my ($arg_href) = @_;

    ## Default(s)
    my $trap_signals_ref;
    my $trap_function;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
        FILEHANDLE       => { required => 1, store => \$FILEHANDLE },
        trap_signals_ref => {
            default     => ["ERR"],
            strict_type => 1,
            store       => \$trap_signals_ref
        },
        trap_function => {
            default     => "error",
            strict_type => 1,
            store       => \$trap_function
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    print $FILEHANDLE "\n## Enable trap for signal(s) "
      . join( " ", @$trap_signals_ref ), "\n";
    print $FILEHANDLE "trap '"
      . $trap_function . "' "
      . join( " ", @$trap_signals_ref ), "\n\n";
}

sub track_progress {

##track_progress

##Function : Output SLURM info on each job via sacct command and write to log file(.status)
##Returns  : ""
##Arguments: $job_ids_ref, $sacct_format_fields_ref, $log_file_ref, $FILEHANDLE
##         : $job_ids_ref             => Job ids
##         : $sacct_format_fields_ref => Format and fields of sacct output
##         : $log_file_ref            => The log file {REF}
##         : $FILEHANDLE              => Sbatch filehandle to write to

    my ($arg_href) = @_;

    ## Default(s)
    my $sacct_format_fields_ref;

    ## Flatten argument(s)
    my $job_ids_ref;
    my $log_file_ref;
    my $FILEHANDLE;

    my $tmpl = {
        job_ids_ref =>
          { default => [], strict_type => 1, store => \$job_ids_ref },
        sacct_format_fields_ref => {
            default => [
                "jobid",     "jobname%50", "account", "partition",
                "alloccpus", "TotalCPU",   "elapsed", "start",
                "end",       "state",      "exitcode"
            ],
            strict_type => 1,
            store       => \$sacct_format_fields_ref
        },
        log_file_ref =>
          { default => \$$, strict_type => 1, store => \$log_file_ref },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    use MIP::Workloadmanager::Slurm qw(slurm_sacct);

    if (@$job_ids_ref) {

        ## Copy array
        my @reformat_sacct_header = @{$sacct_format_fields_ref};

        ## Remove "%digits" from headers
        foreach my $element (@reformat_sacct_header) {

            $element =~ s/%\d+//g;
        }
        my @command = slurm_sacct(
            {
                fields_format_ref => \@{$sacct_format_fields_ref},
                job_ids_ref       => \@{$job_ids_ref},
            }
        );
        print $FILEHANDLE "\t" . join( " ", @command ) . " ";
        print $FILEHANDLE q?| ?;
        print $FILEHANDLE q?perl -nae 'my @headers=(?
          . join( ",", @reformat_sacct_header )
          . q?); if($. == 1) {print "#".join("\t", @headers), "\n"} if ($.>=3 && $F[0]!~/.batch/) {print join("\t", @F), "\n"}' ?;
        print $FILEHANDLE q?> ? . $$log_file_ref . ".status", "\n\n";
    }
}

1;
