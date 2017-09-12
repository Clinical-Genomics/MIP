package MIP::Get::File;

use strict;
use warnings;
use warnings qw{FATAL utf8};
use utf8;    #Allow unicode characters in this script
use open qw{  :encoding(UTF-8) :std};
use charnames qw{ :full :short };
use Carp;
use English qw{-no_match_vars};
use Params::Check qw{check allow last_error};

use Readonly;

BEGIN {
    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{get_file_suffix};
}

## Constants
Readonly my $SPACE => q{ };

sub get_file_suffix {

##get_file_suffix

##Function : Return the current file suffix for this jobid chain or program
##Returns  : "$file_suffix"
##Arguments: $parameter_href, $suffix_key, $job_id_chain, $program_name
##         : $parameter_href => Holds all parameters
##         : $suffix_key     => Suffix key
##         : $job_id_chain   => Job id chain for program
##         : $program_name   => Program name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $suffix_key;
    my $job_id_chain;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        suffix_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$suffix_key
        },
        jobid_chain  => { strict_type => 1, store => \$job_id_chain },
        program_name => { strict_type => 1, store => \$program_name },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $file_suffix;

    ## Jobid chain specific suffix
    if ( defined $job_id_chain ) {

        $file_suffix = $parameter_href->{$suffix_key}{$job_id_chain};
    }
    elsif ( defined $program_name ) {
        ## Program  specific

        $file_suffix = $parameter_href->{$program_name}{$suffix_key};
    }

    ## If suffix was found
    if ( ( defined $file_suffix ) && ($file_suffix) ) {

        return $file_suffix;
    }
    else {
        ## Broadcast no suffic was found
        if ( defined $job_id_chain ) {

            say $log->fatal(
                q{Could not get requested infile_suffix for jobid_chain:}
                  . $job_id_chain );
        }
        elsif ( defined $program_name ) {

            say $log->fatal(
                q{Could not get requested infile_suffix for program:}
                  . $program_name );
        }
        exit 1;
    }
    return;
}

1;
