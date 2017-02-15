package Script::Utils;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

use Exporter qw(import);
 
our @EXPORT_OK = ("help",
		  "set_default_array_parameters",
    );


sub help {

##help

##Function : Print help text and exit with supplied exit code
##Returns  : ""
##Arguments: $USAGE, $exit_code
##         : $USAGE     => Help text
##         : $exit_code => Exit code

    my ($arg_href) = @_;

    ## Default(s)
    my $exit_code;

    ## Flatten argument(s)
    my $USAGE;

    my $tmpl = {
	USAGE => {required => 1, defined => 1, strict_type => 1, store => \$USAGE},
	exit_code => { default => 0,
		       allow => qr/^\d+$/,
		       strict_type => 1, store => \$exit_code},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    print STDOUT $USAGE, "\n";
    exit $exit_code;
}


sub set_default_array_parameters {

##set_default_array_parameters

##Function : Set default for array parameters unless parameter already exists in parameter hash
##Returns  : ""
##Arguments: $parameter_href, $array_parameter_href
##         : $parameter_href       => Parameters hash {REF}
##         : $array_parameter_href => Hold the array parameter defaults as {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $array_parameter_href;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	array_parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$array_parameter_href},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    foreach my $parameter_name (keys %$array_parameter_href) {

	if (! @{ $parameter_href->{$parameter_name} }) {  #Unless parameter already exists

	    $parameter_href->{$parameter_name} = $array_parameter_href->{$parameter_name}{default};
	}
    }
}
