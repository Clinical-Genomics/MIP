package Program::Variantcalling::Peddy;

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
    our @EXPORT_OK = qw(peddy);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

sub peddy {

##peddy

##Function : Perl wrapper for writing peddy recipe to already open $FILEHANDLE or return commands array. Based on peddy 0.2.9.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_prefix_path, $family_file_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $append_stderr_info, $processor_number, $plot
##         : $infile_path         => Infile path
##         : $outfile_prefix_path => Outfile path
##         : $family_file_path    => Family file path
##         : $stderrfile_path     => Stderrfile path
##         : $stdoutfile_path     => Stdoutfile path
##         : $FILEHANDLE          => Filehandle to write to
##         : $append_stderr_info  => Append stderr info to file
##         : $processor_number    => Number of processors to use
##         : $plot                => Generate plots

    my ($arg_href) = @_;

    ## Default(s)
    my $append_stderr_info;
    my $processor_number;
    my $plot;
    
    ## Flatten argument(s)
    my $infile_path;
    my $outfile_prefix_path;
    my $family_file_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = { 
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_prefix_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_prefix_path},
	family_file_path => { required => 1, defined => 1, strict_type => 1, store => \$family_file_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
	processor_number => { default => 4,
			      allow => qr/^\d+$/,
			      strict_type => 1, store => \$processor_number},
	plot => { default => 1,
		  allow => [0, 1],
		  strict_type => 1, store => \$plot},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## peddy
    my @commands = qw(python -m peddy);  #Stores commands depending on input parameters

    ## Options
    if ($processor_number) {

	push(@commands, "--procs ".$processor_number);
    }
    if ($plot) {

	push(@commands, "--plot");
    }
    ## Outfile
    if ($outfile_prefix_path) {

	push(@commands, "--prefix ".$outfile_prefix_path);
    }
    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    if ($family_file_path) {

	push(@commands, $family_file_path);
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
