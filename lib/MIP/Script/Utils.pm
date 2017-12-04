package MIP::Script::Utils;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Cwd;
use File::Spec::Functions qw{ catdir };

## CPANM
use Readonly;
use autodie;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ help set_default_array_parameters create_temp_dir };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq(\n);
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub help {

## Function : Print help text and exit with supplied exit code
## Returns  :
## Arguments: $USAGE     => Help text
##          : $exit_code => Exit code

    my ($arg_href) = @_;

    ## Default(s)
    my $exit_code;

    ## Flatten argument(s)
    my $USAGE;

    my $tmpl = {
        USAGE =>
          { required => 1, defined => 1, strict_type => 1, store => \$USAGE },
        exit_code => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$exit_code,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {*STDOUT} $USAGE;
    exit $exit_code;
}

sub set_default_array_parameters {

## Function : Set default for array parameters unless parameter already exists in parameter hash
## Returns  :
## Arguments: $parameter_href       => Parameters hash {REF}
##          : $array_parameter_href => Hold the array parameter defaults as {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $array_parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        array_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$array_parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter_name ( keys %{$array_parameter_href} ) {

        ## Unless parameter already exists
        if ( not @{ $parameter_href->{$parameter_name} } ) {

            $parameter_href->{$parameter_name} =
              $array_parameter_href->{$parameter_name}{default};
        }
    }
    return;
}

sub create_temp_dir {

## Function : Create a temporary directory and returns the path to it
## Returns  : $temp_dir_path
## Arguments: $directory_path       => Where to create the temporary directory
##          : $directory_base_name  => Base name of directroy
##          : $max_expression_value => Max integrer to add to directory_name
##          : $FILEHANDLE           => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory_path;
    my $directory_base_name;
    my $max_expression_value;
    my $FILEHANDLE;

    my $tmpl = {
        directory_path => {
            defined     => 1,
            default     => cwd(),
            strict_type => 1,
            store       => \$directory_path,
        },
        directory_base_name => {
            defined     => 1,
            default     => q{temp_dir},
            strict_type => 1,
            store       => \$directory_base_name,
        },
        max_expression_value => {
            defined     => 1,
            default     => 1000,
            strict_type => 1,
            store       => \$max_expression_value,
        },
        FILEHANDLE => {
            required => 1,
            defined  => 1,
            store    => \$FILEHANDLE,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_mkdir };

    my $temp_dir_path = catdir( $directory_path,
        $directory_base_name . $UNDERSCORE . int rand $max_expression_value );

    gnu_mkdir(
        {
            indirectory_path => $temp_dir_path,
            parents          => 1,
            FILEHANDLE       => $FILEHANDLE,
        }
    );

    return $temp_dir_path;
}

1;
