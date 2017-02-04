package Script::Utils;

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

use Exporter qw(import);
 
our @EXPORT_OK = ("help",
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

