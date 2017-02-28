package Program::Command::Xargs;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
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
    our @EXPORT_OK = qw(xargs);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

sub xargs {

##xargs

##Function : Perl wrapper for writing xargs view recipe to already open $FILEHANDLE or return command line string. Based on xargs 4.4.2
##Returns  : "$cmd_line"
##Arguments: $FILEHANDLE, $replace_str, $verbose, $max_args, $max_procs, $shell_command, $placeholder_symbol
##         : $FILEHANDLE         => Filehandle to write to
##         : $replace_str        => Replace string.  Enables us to tell xargs where to put the command file lines
##         : $verbose            => Print the command line on the standard error output before executing it
##         : $max_args           => Use at most max-args arguments per command line
##         : $max_procs          => Run up to max-procs processes at a time
##         : $shell_command      => The string following this command will be interpreted as a shell command
##         : $placeholder_symbol => Set placeholder symbol

    my ($arg_href) = @_;

    ## Default(s)
    my $replace_str;
    my $verbose;
    my $max_args;
    my $max_procs;
    my $shell_command;
    my $placeholder_symbol;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = { 
	FILEHANDLE => { store => \$FILEHANDLE},
	replace_str => { default => 0,
			 allow => [0, 1],
			 strict_type => 1, store => \$replace_str},
	verbose => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$verbose},
	max_args => { default => 1,
		      allow => qr/^\d+$/,,
		      strict_type => 1, store => \$max_args},
	max_procs => { default => 1,
		       allow => qr/^\d+$/,
		       strict_type => 1, store => \$max_procs},
	shell_command => { default => 1,
			   strict_type => 1, store => \$shell_command},
	placeholder_symbol => { default => "{}",
				strict_type => 1, store => \$placeholder_symbol},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $cmd_line = "xargs ";
    if($replace_str) {

	$cmd_line .= "-i ";  #replace-str; Enables us to tell xargs where to put the command file lines
    }
    if ($verbose) {

	$cmd_line .= "--verbose ";  #Print the command line on the standard error output before executing it
    }
    if ($max_args) {

	$cmd_line .= "-n ".$max_args." ";  #Use at most max-args arguments per command line
    }
    if ($max_procs) {
	
	$cmd_line .= "-P ".$max_procs." ";  #Run up to max-procs processes at a time
    }
    if ($placeholder_symbol) {

	$cmd_line .= "sh -c".q? "?;  #The string following this command will be interpreted as a shell command
	if ($shell_command) {

	    ##Right trims trailing whitespace from string
	    $shell_command = right_trim({string => $shell_command,
					});
	    #print $FILEHANDLE $shell_command." ";
	    $cmd_line .= $shell_command." ";
	}
	$cmd_line .= $placeholder_symbol.q?" ?;  #Set placeholder and end quotes
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE $cmd_line;
    }
    return $cmd_line;
}


sub right_trim {

##right_trim

##Function : Right trims trailing whitespace from string
##Returns  : ""
##Arguments: $string
##         : $string => The string to trim

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $string;

    my $tmpl = { 
	string => { required => 1, defined => 1, strict_type => 1, store => \$string},
    };
     
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    $string =~ s/\s+$//;  #Remove trailing whitespace
    return $string;
};
