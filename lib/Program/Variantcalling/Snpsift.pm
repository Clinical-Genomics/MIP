package Program::Variantcalling::Snpsift;

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
    our @EXPORT_OK = qw(annotate dbnsfp);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

sub annotate {

##annotate

##Function : Perl wrapper for writing snpsift ann recipe to already open $FILEHANDLE or return commands array. Based on Snpsift 4.2 (build 2015-12-05).
##Returns  : "@commands"
##Arguments: $database_path, $infile_path, $outfile_path, $stderrfile_path, $stdoutfile_path, $config_file_path, $name_prefix, $info, $FILEHANDLE, $append_stderr_info, $verbosity
##         : $database_path        => Database path
##         : $infile_path          => Infile path
##         : $outfile_path         => Outfile path
##         : $stderrfile_path      => Stderrfile path
##         : $stdoutfile_path      => Stdoutfile path
##         : $config_file_path     => Config file path
##         : $name_prefix          => Prepend 'str' to all annotated INFO fields
##         : $info                 => Annotate using a list of info fields (list is a comma separated list of fields)
##         : $FILEHANDLE           => Filehandle to write to
##         : $append_stderr_info   => Append stderr info to file
##         : $verbosity            => Increase output verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $append_stderr_info;
    
    ## Flatten argument(s)
    my $database_path;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $config_file_path;
    my $name_prefix;
    my $info;
    my $FILEHANDLE;
    my $verbosity;

    my $tmpl = { 
	database_path => { required => 1, defined => 1, strict_type => 1, store => \$database_path},
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	config_file_path => { strict_type => 1, store => \$config_file_path},
	name_prefix => { strict_type => 1, store => \$name_prefix},
	info => { strict_type => 1, store => \$info},
	FILEHANDLE => { store => \$FILEHANDLE},
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Snpsift annotate
    my @commands = qw(annotate);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

        push(@commands, "-".$verbosity);
    }
    if ($config_file_path) {

        push(@commands, "-config ".$config_file_path);
    }
    if ($name_prefix) {

        push(@commands, "-name ".$name_prefix);
    }
    if ($info) {

        push(@commands, "-info ".$info);
    }
    if ($database_path) {

        push(@commands, $database_path);
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


sub dbnsfp {

##dbnsfp

##Function : Perl wrapper for writing snpsift dbnsfp recipe to already open $FILEHANDLE or return commands array. Based on Snpsift 4.2 (build 2015-12-05).
##Returns  : "@commands"
##Arguments: $annotate_fields_ref, $database_path, $infile_path, $outfile_path, $stderrfile_path, $stdoutfile_path, $config_file_path, $FILEHANDLE, $append_stderr_info, $verbosity
##         : annotate_fields_ref   => Add annotations for list of fields
##         : $database_path        => Database path
##         : $infile_path          => Infile path
##         : $outfile_path         => Outfile path
##         : $stderrfile_path      => Stderrfile path
##         : $stdoutfile_path      => Stdoutfile path
##         : $config_file_path     => Config file path
##         : $FILEHANDLE           => Filehandle to write to
##         : $append_stderr_info   => Append stderr info to file
##         : $verbosity            => Increase output verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $append_stderr_info;
    
    ## Flatten argument(s)
    my $annotate_fields_ref;
    my $database_path;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $config_file_path;
    my $FILEHANDLE;
    my $verbosity;

    my $tmpl = {
	annotate_fields_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$annotate_fields_ref},
	database_path => { required => 1, defined => 1, strict_type => 1, store => \$database_path},
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

    ## Snpsift dbnsfp
    my @commands = qw(dbnsfp);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

        push(@commands, "-".$verbosity);
    }
    if ($config_file_path) {

        push(@commands, "-config ".$config_file_path);
    }
    if ($database_path) {

        push(@commands, "-db ".$database_path);
    }
    if (@$annotate_fields_ref) {

        push(@commands, "-f ".join(",", @$annotate_fields_ref));
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
