package Program::Command::Gnu;

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
    our @EXPORT_OK = qw(cp rm split xargs);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub cp {

##cp

##Function : Perl wrapper for writing cp recipe to already open $FILEHANDLE or return commands array. Based on cp 8.4
##Returns  : "@commands"
##Arguments: $preserve_attributes_ref, $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $preserve, $recursive, $force, $verbose
##         : $preserve_attributes_ref => Preserve the specified attributes (default:mode,ownership,timestamps), if possible additional attributes: context, links, xattr, all
##         : $FILEHANDLE              => Filehandle to write to
##         : $infile_path             => Infile path
##         : $outfile_path            => Outfile path
##         : $stderrfile_path         => Stderrfile path
##         : $preserve                => Same as --preserve=mode,ownership,timestamps
##         : $recursive               => Copy directories recursively
##         : $force                   => If an existing destination file cannot be opened, remove it and try again
##         : $verbose                 => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $preserve;
    my $recursive;
    my $force;
    my $verbose;

    ## Flatten argument(s)
    my $preserve_attributes_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;

    my $tmpl = { 
	preserve_attributes_ref => { default => [], strict_type => 1, store => \$preserve_attributes_ref},
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	recursive => { default => 0,
		       strict_type => 1, store => \$recursive},
	force => { default => 0,
		   strict_type => 1, store => \$force},
	preserve => { default => 0,
		      strict_type => 1, store => \$preserve},
	verbose => { default => 0,
		     strict_type => 1, store => \$verbose},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cp
    my @commands = qw(cp);  #Stores commands depending on input parameters

    if(@$preserve_attributes_ref) {

	push(@commands, "--preserve=".join(",", @$preserve_attributes_ref));  #Preserve the specified attributes
    }
    elsif($preserve) {

	push(@commands, "-p");
    }
    if ($recursive) {

	push(@commands, "--recursive");
    }
    if ($force) {

	push(@commands, "--force");
    }
    if ($verbose) {

	push(@commands, "--verbose");  #Explain what is being done
    }
    push(@commands, $infile_path);
    push(@commands, $outfile_path);

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub rm {

##rm

##Function : Perl wrapper for writing rm recipe to already open $FILEHANDLE or return commands array. Based on rm 8.4
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $stderrfile_path, $recursive, $force, $verbose
##         : $preserve_attributes_ref => Preserve the specified attributes (default:mode,ownership,timestamps), if possible additional attributes: context, links, xattr, all
##         : $FILEHANDLE              => Filehandle to write to
##         : $infile_path             => Infile path
##         : $stderrfile_path         => Stderrfile path
##         : $recursive               => Copy directories recursively
##         : $force                   => If an existing destination file cannot be opened, remove it and try again
##         : $verbose                 => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $recursive;
    my $force;
    my $verbose;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $stderrfile_path;

    my $tmpl = { 
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	recursive => { default => 0,
		       strict_type => 1, store => \$recursive},
	force => { default => 0,
		   strict_type => 1, store => \$force},
	verbose => { default => 0,
		     strict_type => 1, store => \$verbose},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## rm
    my @commands = qw(rm);  #Stores commands depending on input parameters

    if ($recursive) {

	push(@commands, "--recursive");
    }
    if ($force) {

	push(@commands, "--force");
    }
    if ($verbose) {

	push(@commands, "--verbose");  #Explain what is being done
    }

    ## Infile
    push(@commands, $infile_path);

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub split {

##split

##Function : Perl wrapper for writing split recipe to $FILEHANDLE or return commands array. Based on split 8.4.
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $prefix, $lines, $suffix_length, $numeric_suffixes, $quiet, $verbose
##         : $FILEHANDLE       => Filehandle to write to
##         : $infile_path      => Infile path
##         : $prefix           => Prefix of output files
##         : $lines            => Put number lines per output file
##         : $suffix_length    => Use suffixes of length N
##         : $numeric_suffixes => Use numeric suffixes instead of alphabetic
##         : $quiet            => Suppress all warnings
##         : $verbose          => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $quiet;
    my $verbose;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $prefix;
    my $lines;
    my $suffix_length;
    my $numeric_suffixes;

    my $tmpl = {
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	prefix => { strict_type => 1, store => \$prefix},
	lines => { allow => qr/^\d+$/,
		   strict_type => 1, store => \$lines},
	suffix_length => { allow => qr/^\d+$/,
			   strict_type => 1, store => \$suffix_length},
	numeric_suffixes => { default => 0,
			      allow => [0, 1],
			      strict_type => 1, store => \$numeric_suffixes},
	quiet => { default => 0,
		   allow => [0, 1],
		   strict_type => 1, store => \$quiet},
	verbose => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$verbose},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Split
    my @commands = qw(split);  #Stores commands depending on input parameters

    ## Options
    if ($lines) {

	push(@commands, "--lines=".$lines);
    }
    if ($numeric_suffixes) {

	push(@commands, "--numeric-suffixes");
    }
    if ($numeric_suffixes) {
	
	push(@commands, "--suffix-length=".$numeric_suffixes);
    }
    if ($quiet) {

	push(@commands, "--quiet");
    }
    if ($verbose) {

	push(@commands, "--verbose");
    }

    ## FILE
    push(@commands, $infile_path);

    if ($prefix) {

	push(@commands, $prefix);
    }

    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub xargs {

##xargs

##Function : Perl wrapper for writing xargs recipe to already open $FILEHANDLE or return command line string. Based on xargs 4.4.2
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


1;
