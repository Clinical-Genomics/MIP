package File::Parse::Parse;

use warnings qw( FATAL utf8 );
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case
use Cwd qw(abs_path);  #Import absolute path function

use Exporter qw(import);
 
our @EXPORT_OK = ("find_absolute_path",
    );


sub find_absolute_path {

##find_absolute_path

##Function : Find aboslute path for supplied path or croaks and exists if path does not exists
##Returns  : "$path - absolute path"
##Arguments: $path, $parameter_name
##         : $path           => Supplied path to be updated/evaluated
##         : $parameter_name => Parameter to be evaluated
##         : $log            => Log object to write to if supplied

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $path;
    my $parameter_name;
    my $log;

    my $tmpl = {
	path => { required => 1, defined => 1, store => \$path},
	parameter_name => { required => 1, defined => 1, strict_type => 1, store => \$parameter_name},
	log => { store => \$log},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $original_path = $path;

    $path = abs_path($path);

    unless(defined($path)) {

	if ( (defined($log)) && ($log) ) {

	    $log->fatal("Could not find absolute path for ".$parameter_name.": ".$original_path.". Please check the supplied path!\n");
	    exit 1;
	}
	else {
	    
	    warn("Could not find absolute path for ".$parameter_name.": ".$original_path.". Please check the supplied path!\n");
	    exit;
	}
    }
    return $path;
}
