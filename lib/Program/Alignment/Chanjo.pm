package Program::Alignment::Chanjo;

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
    our @EXPORT_OK = qw(sex);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub sex {

##sex

##Function : Perl wrapper for writing chanjo sex recipe to $FILEHANDLE. Based on chanjo 4.0.0
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $log_file_path, $chr_prefix, $log_level
##         : $FILEHANDLE      => Sbatch filehandle to write to
##         : $infile_path     => Infile path
##         : $outfile_path    => Outfile path
##         : $stderrfile_path => Stderrfile path
##         : $log_file_path   => Log file path
##         : $chr_prefix      => Chromosome prefix
##         : $log_level       => Level of logging

    my ($arg_href) = @_;

    ## Default(s)
    my $log_level;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $log_file_path;
    my $chr_prefix;
    
    my $tmpl = {
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	log_file_path => { strict_type => 1, store => \$log_file_path },
	chr_prefix => { allow => [undef, "chr"],
			strict_type => 1, store => \$chr_prefix },
	log_level => { default => "INFO",
		       allow => ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
		       strict_type => 1, store => \$log_level },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Chanjo
    my @commands = qw(chanjo);  #Stores commands depending on input parameters

    ## Chanjo main options
    if ($log_level) {

	push(@commands, "--log-level ".$log_level);
    }
    if ($log_file_path) {

	push(@commands, "--log-file ".$log_file_path);
    }

    push(@commands, "sex");

    ## Options
    if ($chr_prefix) {

	push(@commands, "--prefix ".$chr_prefix);
    }
    ##Infile
    push(@commands, $infile_path);

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


1;
