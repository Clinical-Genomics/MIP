package MIP::Check::Parameter;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };

# Allow unicode characters in this script
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{ check allow last_error };

use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ check_allowed_array_values check_allowed_temp_directory};
}

##Constants
Readonly my $NEWLINE => qq{\n};

sub check_allowed_array_values {

##check_allowed_array_values

##Function : Check that the array values are allowed
##Returns  : ""
##Arguments: $allowed_values_ref, $values_ref
##         : $allowed_values_ref => Allowed values for parameter {REF}
##         : $values_ref         => Values for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $allowed_values_ref;
    my $values_ref;

    my $tmpl = {
        allowed_values_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$allowed_values_ref
        },
        values_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$values_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %is_allowed;

    # Remake allowed values into keys in is_allowed hash
    map { $is_allowed{$_} = undef } @{$allowed_values_ref};

  VALUES:
    foreach my $value ( @{$values_ref} ) {

        # Test if value is allowed
        if ( not exists $is_allowed{$value} ) {

            return 0;
        }
    }

    # All ok
    return 1;
}

sub check_allowed_temp_directory {

##check_allowed_temp_directory

##Function : Check that the temp directory value is allowed
##Returns  : ""
##Arguments: $temp_directory
##         : $temp_directory => Temp directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $temp_directory;

    my $tmpl = {
        temp_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$temp_directory
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %is_not_allowed = (
        q{/scratch}  => undef,
        q{/scratch/} => undef,
    );

    # Test if value is allowed
    if ( exists $is_not_allowed{$temp_directory} ) {

        return 0;
    }

    # All ok
    return 1;
}

1;
