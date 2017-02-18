package File::Format::Shell;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Inherit from Exporter to export functions and variables
    our @ISA = qw(Exporter);

    # Functions and variables which are exported by default
    our @EXPORT = qw();

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = ("create_bash_file",
		      "create_housekeeping_function",
		      "create_trap_function",
		      "enable_trap",
		      "clear_trap",
		      "track_progress",
	);
}

use Cwd;
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catfile catdir devnull);
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub create_bash_file {

##create_bash_file

##Function : Create bash file with header
##Returns  : "$FILEHANDLE"
##Arguments: $file_name, $directory_remove, $log, $trap_signals_ref, $trap_function
##         : $file_name        => File name
##         : $directory_remove => Directory to remove when caught by trap function
##         : $log              => Log object to write to
##         : $trap_signals_ref => Array with signals to clear trap for {REF}
##         : $trap_function    => Trap function argument

    my ($arg_href) = @_;

    ## Default(s)
    my $trap_signals_ref;
    my $trap_function;

    ## Flatten argument(s)
    my $file_name;
    my $directory_remove;
    my $log;    

    my $tmpl = {
	file_name => { required => 1, defined => 1, strict_type => 1, store => \$file_name},
	log => { store => \$log},
	directory_remove => { allow => qr/^\.\S+$/,
			      strict_type => 1, store => \$directory_remove},
	trap_signals_ref => { default => ["ERR"],
			      strict_type => 1, store => \$trap_signals_ref},
	trap_function => { default => "error",
			   strict_type => 1, store => \$trap_function},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $pwd = cwd();

    ## Open batch file with supplied log
    if ( (defined($log)) && ($log) ) {
	
	open ($FILEHANDLE, ">", catfile($pwd, $file_name)) or $log->logdie("Cannot write to '".catfile($pwd, $file_name)."' :".$!."\n");
    }
    else {

	open ($FILEHANDLE, ">", catfile($pwd, $file_name)) or die("Cannot write to '".catfile($pwd, $file_name)."' :".$!."\n");
    }

    say $FILEHANDLE "#!".catfile( dirname( dirname( devnull() ) ) ).catfile("usr", "bin", "env", "bash"), "\n";

    ## Create housekeeping function which removes entire directory when finished
    create_housekeeping_function({directory_remove => $directory_remove,
				  FILEHANDLE => $FILEHANDLE,
				 });

    ## Create error handling function and trap
    create_trap_function({FILEHANDLE => $FILEHANDLE,
			 });

    if ( (defined($log)) && ($log) ) {

	$log->info("Created bash file: '".catfile($pwd, $file_name), "'\n");
    }
    else {

	say STDERR "Created bash file: '".catfile($pwd, $file_name), "'";
    }
    return $FILEHANDLE;
}


sub create_housekeeping_function {

##create_housekeeping_function

##Function : Create housekeeping function which removes entire directory when finished
##Returns  : ""
##Arguments: $job_id_href, $log_file_ref, $FILEHANDLE, $directory_remove, $trap_signals_ref, $trap_function
##         : $job_id_href      => Job_id hash {REF}
##         : $log_file_ref     => Log file to write job_id progress to {REF}
##         : $FILEHANDLE       => Filehandle to write to
##         : $directory_remove => Directory to remove when caught by trap function
##         : $trap_signals_ref => Array with signals to enable trap for {REF}
##         : $trap_function    => The trap function argument

    my ($arg_href) = @_;

    ## Default(s)
    my $trap_signals_ref;
    my $trap_function;

    ## Flatten argument(s)
    my $job_id_href;
    my $log_file_ref;
    my $FILEHANDLE;
    my $directory_remove;

    my $tmpl = {
	job_id_href => { default => {}, strict_type => 1, store => \$job_id_href},
	log_file_ref => { default => \$$, strict_type => 1, store => \$log_file_ref},
	FILEHANDLE => { required => 1, store => \$FILEHANDLE},
	directory_remove => {strict_type => 1, store => \$directory_remove},
	trap_signals_ref => { default => ["EXIT", "TERM", "INT"],
			      strict_type => 1, store => \$trap_signals_ref},
	trap_function => { default => "finish",
			   strict_type => 1, store => \$trap_function},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Create housekeeping function and trap
    say $FILEHANDLE q?finish() {?, "\n";

    if ( (defined($directory_remove)) && ($directory_remove) ) {
	
	say $FILEHANDLE "\t".q?## Perform exit housekeeping?;
	say $FILEHANDLE "\t".q?rm -rf ?.$directory_remove;
    }
    if ( (defined($job_id_href)) && (keys %$job_id_href)
	&& (defined($$log_file_ref)) && ($$log_file_ref) ) {

	## Output SLURM info on each job via sacct command and write to log file(.status)
	track_progress({job_id_href => $job_id_href,
			FILEHANDLE => $FILEHANDLE,
			log_file_ref => $log_file_ref,
		       });
    }
    say $FILEHANDLE q?}?;

    ## Enable trap function with trap signal(s)
    enable_trap({FILEHANDLE => $FILEHANDLE,
		 trap_signals_ref => \@{ $trap_signals_ref },
		 trap_function => $trap_function,
		});
}


sub create_trap_function {

##create_trap_function

##Function : Create error handling function and trap
##Returns  : ""
##Arguments: $job_id_href, $log_file_ref, $FILEHANDLE, $trap_signals_ref, $trap_function
##         : $job_id_href      => Job_id hash {REF}
##         : $log_file_ref     => Log file to write job_id progress to {REF}
##         : $FILEHANDLE       => Filehandle to write to
##         : $trap_signals_ref => Array with signals to enable trap for {REF}
##         : $trap_function    => The trap function argument

    my ($arg_href) = @_;

    ## Default(s)
    my $trap_signals_ref;
    my $trap_function;

    ## Flatten argument(s)
    my $job_id_href;
    my $log_file_ref;
    my $FILEHANDLE;

    my $tmpl = {
	job_id_href => { default => {}, strict_type => 1, store => \$job_id_href},
	log_file_ref => { default => \$$, strict_type => 1, store => \$log_file_ref},
	FILEHANDLE => { required => 1, store => \$FILEHANDLE},
	trap_signals_ref => { default => ["ERR"],
			      strict_type => 1, store => \$trap_signals_ref},
	trap_function => { default => "error",
			   strict_type => 1, store => \$trap_function},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Create error handling function and trap
    say $FILEHANDLE $trap_function.q?() {?, "\n";
    say $FILEHANDLE "\t".q?## Display error message and exit?;
    say $FILEHANDLE "\t".q{ret="$?"};
    say $FILEHANDLE "\t".q?echo "${PROGNAME}: ${1:-"Unknown Error - ExitCode="$ret}" 1>&2?, "\n";
    say $FILEHANDLE "\t".q?exit 1?;

    if ( (defined($job_id_href)) && (keys %$job_id_href)
	&& (defined($$log_file_ref)) && ($$log_file_ref) ) {

	## Output SLURM info on each job via sacct command and write to log file(.status)
	track_progress({job_id_href => $job_id_href,
			FILEHANDLE => $FILEHANDLE,
			log_file_ref => $log_file_ref,
		       });
    }
    say $FILEHANDLE q?}?;

    ## Enable trap function with trap signal(s)
    enable_trap({FILEHANDLE => $FILEHANDLE,
		 trap_signals_ref => \@{ $trap_signals_ref },
		 trap_function => $trap_function,
		});
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
	FILEHANDLE => { required => 1, store => \$FILEHANDLE},
	trap_signals_ref => { default => ["ERR"],
			      strict_type => 1, store => \$trap_signals_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Clear trap for signal ERR
    say $FILEHANDLE "\n## Clear trap for signal(s) ".join(" ", @$trap_signals_ref);
    say $FILEHANDLE "trap - ".join(" ", @$trap_signals_ref);
    say $FILEHANDLE "trap", "\n";
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
	FILEHANDLE => { required => 1, store => \$FILEHANDLE},
	trap_signals_ref => { default => ["ERR"],
			      strict_type => 1, store => \$trap_signals_ref},
	trap_function => { default => "error",
			   strict_type => 1, store => \$trap_function},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    say $FILEHANDLE "\n## Enable cleared trap for signal(s) ".join(" ", @$trap_signals_ref)." ";
    say $FILEHANDLE "trap ".$trap_function." ".join(" ", @$trap_signals_ref), "\n";
}


sub track_progress {

##track_progress

##Function : Output SLURM info on each job via sacct command and write to log file(.status)
##Returns  : ""
##Arguments: $job_id_href, $log_file_ref, $FILEHANDLE
##         : $job_id_href  => The job_id hash {REF}
##         : $log_file_ref => The log file {REF}
##         : $FILEHANDLE   => Sbatch filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_href;
    my $log_file_ref;
    my $FILEHANDLE;

    my $tmpl = {
	job_id_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$job_id_href},
	log_file_ref => { default => \$$, strict_type => 1, store => \$log_file_ref},
	FILEHANDLE => { store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    if (keys %$job_id_href) {

	print $FILEHANDLE "\t".q?sacct --format=jobid,jobname%50,account,partition,alloccpus,TotalCPU,elapsed,start,end,state,exitcode -j ?;
	print $FILEHANDLE join(',', @{ $job_id_href->{PAN}{PAN} }), " ";
	print $FILEHANDLE q?| ?;
	print $FILEHANDLE q?perl -nae 'my @headers=(jobid,jobname,account,partition,alloccpus,TotalCPU,elapsed,start,end,state,exitcode); if($. == 1) {print "#".join("\t", @headers), "\n"} if ($.>=3 && $F[0]!~/.batch/) {print join("\t", @F), "\n"}' ?;
	say $FILEHANDLE q?> ?.$$log_file_ref.".status";
    }
}


1;
