package Program::Qc::Mip;

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
    our @EXPORT_OK = qw(qccollect);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

sub qccollect {

##qccollect

##Function : Perl wrapper for writing qccollect recipe to already open $FILEHANDLE or return commands array. Based on qccollect 2.0.2.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $regexp_file_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $append_stderr_info, $skip_evaluation
##         : $infile_path        => Infile path
##         : $outfile_path       => Outfile path
##         : $regexp_file_path   => Regular expression file
##         : $stderrfile_path    => Stderrfile path
##         : $stdoutfile_path    => Stdoutfile path
##         : $FILEHANDLE         => Filehandle to write to
##         : $append_stderr_info => Append stderr info to file
##         : $skip_evaluation    => Skip evaluation step

    my ($arg_href) = @_;

    ## Default(s)
    my $append_stderr_info;
    my $skip_evaluation;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $regexp_file_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = { 
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path},
	regexp_file_path => { required => 1, defined => 1, strict_type => 1, store => \$regexp_file_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
	skip_evaluation => { default => 0,
			     allow => [0, 1],
			     strict_type => 1, store => \$skip_evaluation},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Qccollect
    my @commands = qw(qccollect);  #Stores commands depending on input parameters

    ## Options
    if ($regexp_file_path) {

	push(@commands, "--regexp_file ".$regexp_file_path);
    }
    if ($skip_evaluation) {

	push(@commands, "--skip_evaluation ".$skip_evaluation);
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    ## Outfile
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);
    }
    if ($stdoutfile_path) {

	push(@commands, "1> ".$stdoutfile_path);  #Redirect stdout to program specific stdout file
    }
    if ($stderrfile_path) {

	if ($append_stderr_info) {

	    push(@commands, "2>> ".$stderrfile_path);  #Redirect and append stderr output to program specific stderr file
	}
	else {
	    
	    push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
	}
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
