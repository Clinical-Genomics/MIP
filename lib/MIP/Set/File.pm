package MIP::Set::File;

use strict;
use warnings;
use warnings qw{FATAL utf8};
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

BEGIN {

    use base qw{Exporter};
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{set_file_suffix};
}

sub set_file_suffix {

##set_file_suffix

##Function : Set the current file suffix for this job id chain
##Returns  : "$file_suffix"
##Arguments: $parameter_href, $suffix_key, $job_id_chain, $file_suffix
##         : $parameter_href => Holds all parameters
##         : $suffix_key     => Suffix key
##         : $job_id_chain   => Job id chain for program
##         : $file_suffix    => File suffix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $suffix_key;
    my $job_id_chain;
    my $file_suffix;

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
        job_id_chain => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$job_id_chain
        },
        file_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_suffix
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $parameter_href->{$suffix_key}{$job_id_chain} = $file_suffix;

    return $file_suffix;
}

1;
