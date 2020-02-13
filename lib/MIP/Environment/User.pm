package MIP::Environment::User;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_email_address };
}

sub check_email_address {

## Function : Check the syntax of the email adress is valid and has a mail host.
## Returns  :
## Arguments: $email => Email adress

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $email;

    my $tmpl = {
        email => {
            required    => 1,
            store       => \$email,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Email::Valid;

    return if ( not defined $email );

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check syntax and mail host
    my $address = Email::Valid->address(
        -address => $email,
        -mxcheck => 1,
    );
    if ( not defined $address ) {

        $log->fatal( q{The supplied email: }
              . $email
              . q{ seem to be malformed according to }
              . Email::Valid->details() );
        exit 1;
    }
    return 1;
}

1;
