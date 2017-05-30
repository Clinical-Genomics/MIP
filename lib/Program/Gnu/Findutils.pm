package Program::Gnu::Findutils;

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

##Function : Perl wrapper for writing xargs recipe to already open $FILEHANDLE or return command line string. Based on xargs 4.4.2
##Returns  : "@commands"
##Arguments: $shell_commands_ref, $FILEHANDLE, $replace_str, $verbose, $max_args, $max_procs, $placeholder_symbol
##         : $shell_commands_ref => The string following this command will be interpreted as a shell command {REF}
##         : $FILEHANDLE         => Filehandle to write to
##         : $replace_str        => Replace string.  Enables us to tell xargs where to put the command file lines
##         : $verbose            => Print the command line on the standard error output before executing it
##         : $max_args           => Use at most max-args arguments per command line
##         : $max_procs          => Run up to max-procs processes at a time
##         : $placeholder_symbol => Set placeholder symbol

    my ($arg_href) = @_;

    ## Default(s)
    my $replace_str;
    my $verbose;
    my $max_args;
    my $max_procs;
    my $placeholder_symbol;

    ## Flatten argument(s)
    my $shell_commands_ref;
    my $FILEHANDLE;

    my $tmpl = { 
	shell_commands_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$shell_commands_ref},
	FILEHANDLE => { store => \$FILEHANDLE },
	replace_str => { default => 0,
			 allow => [0, 1],
			 strict_type => 1, store => \$replace_str },
	verbose => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$verbose },
	max_args => { default => 1,
		      allow => qr/^\d+$/,,
		      strict_type => 1, store => \$max_args },
	max_procs => { default => 1,
		       allow => qr/^\d+$/,
		       strict_type => 1, store => \$max_procs },
	placeholder_symbol => { default => "{}",
				strict_type => 1, store => \$placeholder_symbol },
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Xargs
    my @commands = qw(xargs);  #Stores commands depending on input parameters

    if($replace_str) {

	push(@commands, "-i");  #replace-str; Enables us to tell xargs where to put the command file lines
    }
    if ($verbose) {

	push(@commands, "--verbose");  #Print the command line on the standard error output before executing it
    }
    if ($max_args) {

	push(@commands, "-n ".$max_args);  #Use at most max-args arguments per command line
    }
    if ($max_procs) {
	
	push(@commands, "-P ".$max_procs);  #Run up to max-procs processes at a time
    }
    if ($placeholder_symbol) {

	push(@commands, "sh -c".q? "?);  #The string following this command will be interpreted as a shell command

	if ($shell_commands_ref) {

	    push(@commands, @$shell_commands_ref);
	}
	push(@commands, $placeholder_symbol.q?" ?);  #Set placeholder and end quotes
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
