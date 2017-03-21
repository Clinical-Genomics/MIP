package Program::Gnu::Coreutils;

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
    our @EXPORT_OK = qw(cp rm mv mkdir cat echo split);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub cp {

##cp

##Function : Perl wrapper for writing cp recipe to already open $FILEHANDLE or return commands array. Based on cp 8.4
##Returns  : "@commands"
##Arguments: $preserve_attributes_ref, $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $preserve, $recursive, $force, $verbose
##         : $preserve_attributes_ref => Preserve the specified attributes (default:mode,ownership,timestamps), if possible additional attributes: context, links, xattr, all
##         : $infile_path             => Infile path
##         : $outfile_path            => Outfile path
##         : $stderrfile_path         => Stderrfile path
##         : $FILEHANDLE              => Filehandle to write to
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
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = { 
	preserve_attributes_ref => { default => [], strict_type => 1, store => \$preserve_attributes_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	recursive => { default => 0,
		       allow => [0, 1],
		       strict_type => 1, store => \$recursive },
	force => { default => 0,
		   allow => [0, 1],
		   strict_type => 1, store => \$force },
	preserve => { default => 0,
		      allow => [0, 1],
		      strict_type => 1, store => \$preserve },
	verbose => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$verbose },
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


sub mv {

##mv

##Function : Perl wrapper for writing mv recipe to already open $FILEHANDLE or return commands array. Based on mv 8.4
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $force, $verbose
##         : $infile_path             => Infile path
##         : $outfile_path            => Outfile path
##         : $stderrfile_path         => Stderrfile path
##         : $FILEHANDLE              => Filehandle to write to
##         : $force                   => If an existing destination file cannot be opened, remove it and try again
##         : $verbose                 => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $force;
    my $verbose;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = { 
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	force => { default => 0,
		   allow => [0, 1],
		   strict_type => 1, store => \$force },
	verbose => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$verbose },
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## mv
    my @commands = qw(mv);  #Stores commands depending on input parameters

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
##Arguments: $infile_path, $stderrfile_path, $FILEHANDLE, $recursive, $force, $verbose
##         : $infile_path     => Infile path
##         : $stderrfile_path => Stderrfile path
##         : $FILEHANDLE      => Filehandle to write to
##         : $recursive       => Copy directories recursively
##         : $force           => If an existing destination file cannot be opened, remove it and try again
##         : $verbose         => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $recursive;
    my $force;
    my $verbose;

    ## Flatten argument(s)
    my $infile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = { 
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	recursive => { default => 0,
		       allow => [0, 1],
		       strict_type => 1, store => \$recursive },
	force => { default => 0,
		   allow => [0, 1],
		   strict_type => 1, store => \$force },
	verbose => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$verbose },
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


sub mkdir {

##mkdir

##Function : Perl wrapper for writing mkdir recipe to already open $FILEHANDLE or return commands array. Based on mkdir 8.4
##Returns  : "@commands"
##Arguments: $indirectory_path, $stderrfile_path, $FILEHANDLE, $parents, $verbose
##         : $indirectory_path => Infile path
##         : $stderrfile_path  => Stderrfile path
##         : $FILEHANDLE       => Filehandle to write to
##         : $parents          => No error if existing, make parent directories as needed
##         : $verbose          => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $verbose;
    my $parents;

    ## Flatten argument(s)
    my $indirectory_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = { 
	indirectory_path => { required => 1, defined => 1, strict_type => 1, store => \$indirectory_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	parents => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$parents },
	verbose => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$verbose },
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## mkdir
    my @commands = qw(mkdir);  #Stores commands depending on input parameters

    if ($parents) {

	push(@commands, "--parents");  #Make parent directories as needed
    }
    if ($verbose) {

	push(@commands, "--verbose");  #Explain what is being done
    }

    ## Indirectory
    push(@commands, $indirectory_path);

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub cat {

##cat

##Function : Perl wrapper for writing cat recipe to already open $FILEHANDLE or return commands array. Based on cat 8.4
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $outfile_path, $stderrfile_path, $FILEHANDLE
##         : $infile_paths_ref => Infile paths {REF}
##         : $outfile_path     => Outfile path
##         : $stderrfile_path  => Stderrfile path
##         : $FILEHANDLE       => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = { 
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cat
    my @commands = qw(cat);  #Stores commands depending on input parameters

    ## Infiles
    push(@commands, join(" ", @$infile_paths_ref));

    ## Outfile
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub echo {

##echo

##Function : Perl wrapper for writing echo recipe to already open $FILEHANDLE or return commands array. Based on echo 8.4
##Returns  : "@commands"
##Arguments: $strings_ref, $outfile_path, $stderrfile_path, $FILEHANDLE, $enable_interpretation, $no_trailing_newline
##         : $strings_ref           => Strings to echo {REF}
##         : $outfile_path          => Outfile path
##         : $stderrfile_path       => Stderrfile path
##         : $FILEHANDLE            => Filehandle to write to
##         : $enable_interpretation => Enable interpretation of backslash escapes
##         : $no_trailing_newline   => Do not output the trailing newline

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $strings_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $enable_interpretation;
    my $no_trailing_newline;

    my $tmpl = { 
	strings_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$strings_ref },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	enable_interpretation => { default => 0,
				   allow => [0, 1],
				   strict_type => 1, store => \$enable_interpretation },
	no_trailing_newline => { default => 0,
				   allow => [0, 1],
				   strict_type => 1, store => \$no_trailing_newline },
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Echo
    my @commands = qw(echo);  #Stores commands depending on input parameters

    ##Options
    if ($enable_interpretation) {

	push(@commands, "-e");
    }
    if ($no_trailing_newline) {

	push(@commands, "-n");
    }

    ## Strings
    push(@commands, q?" ?.join(" ", @$strings_ref).q? "?);

    ## Outfile
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);
    }
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
##Arguments: $infile_path, $FILEHANDLE, $prefix, $lines, $suffix_length, $numeric_suffixes, $quiet, $verbose
##         : $infile_path      => Infile path
##         : $FILEHANDLE       => Filehandle to write to
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
    my $infile_path;
    my $FILEHANDLE;
    my $prefix;
    my $lines;
    my $suffix_length;
    my $numeric_suffixes;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	prefix => { strict_type => 1, store => \$prefix },
	lines => { allow => qr/^\d+$/,
		   strict_type => 1, store => \$lines },
	suffix_length => { allow => qr/^\d+$/,
			   strict_type => 1, store => \$suffix_length },
	numeric_suffixes => { default => 0,
			      allow => [0, 1],
			      strict_type => 1, store => \$numeric_suffixes },
	quiet => { default => 0,
		   allow => [0, 1],
		   strict_type => 1, store => \$quiet },
	verbose => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$verbose },
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

    ## Infile
    push(@commands, $infile_path);

    if ($prefix) {

	push(@commands, $prefix);
    }

    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
