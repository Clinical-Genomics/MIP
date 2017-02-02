package Check::Check_modules;

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Params::Check qw[check allow last_error];
 
use Exporter qw(import);
 
our @EXPORT_OK = ("check_modules",
    );

sub check_modules {

##check_modules

##Function : Evaluate that all modules required are installed
##Returns  : ""
##Arguments: $modules_ref, $program_name
##         : $modules_ref  => Array of module names
##         : $program_name => Program name 

    local $Params::Check::PRESERVE_CASE = 1;

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $modules_ref;
    my $program_name;

    my $tmpl = {
	modules_ref => { required => 1, default => [], strict_type => 1, store => \$modules_ref},
	program_name => { required => 1, defined => 1, strict_type => 1, store => \$program_name},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    foreach my $module (@$modules_ref) {

	$module =~s/::/\//g;  #Replace "::" with "/" since the automatic replacement magic only occurs for barewords.
	$module .= ".pm";  #Add perl module ending for the same reason

	eval {

	    require $module;
	};
	if($@) {

	    warn("NOTE: ".$module." not installed - Please install to run ".$program_name.".\n");
	    warn("NOTE: Aborting!\n");
	    exit 1;
	}
    }
}


return 1;
