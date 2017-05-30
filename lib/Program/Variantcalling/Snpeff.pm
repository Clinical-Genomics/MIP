package Program::Variantcalling::Snpeff;

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
    our @EXPORT_OK = qw(ann);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

sub ann {

##ann

##Function : Perl wrapper for writing snpeff ann recipe to already open $FILEHANDLE or return commands array. Based on SnpEff 4.2 (build 2015-12-05).
##Returns  : "@commands"
##Arguments: $genome_build_version, $infile_path, $outfile_path, $stderrfile_path, $stdoutfile_path, $config_file_path, $FILEHANDLE, $append_stderr_info, $verbosity
##         : $infile_path          => Infile path
##         : $outfile_path         => Outfile path
##         : $stderrfile_path      => Stderrfile path
##         : $stdoutfile_path      => Stdoutfile path
##         : $genome_build_version => Genome build version
##         : $config_file_path     => Config file path
##         : $FILEHANDLE           => Filehandle to write to
##         : $append_stderr_info   => Append stderr info to file
##         : $verbosity            => Increase output verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $append_stderr_info;
    
    ## Flatten argument(s)
    my $genome_build_version;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $config_file_path;
    my $FILEHANDLE;
    my $verbosity;

    my $tmpl = { 
	genome_build_version => { required => 1, defined => 1, strict_type => 1, store => \$genome_build_version},
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	config_file_path => { strict_type => 1, store => \$config_file_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Snpeff ann
    my @commands = qw(ann);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

        push(@commands, "-".$verbosity);
    }
    if ($genome_build_version) {

        push(@commands, $genome_build_version);
    }
    if ($config_file_path) {

        push(@commands, "-config ".$config_file_path);
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
