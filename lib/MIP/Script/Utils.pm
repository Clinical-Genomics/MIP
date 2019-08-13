package MIP::Script::Utils;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPANM
use Readonly;
use autodie;

## MIPs lib/
use MIP::Constants qw{ $COLON $DOT $NEWLINE $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      create_temp_dir
      help
      print_parameter_defaults };
}

sub help {

## Function : Print help text and exit with supplied exit code
## Returns  :
## Arguments: $exit_code => Exit code
##          : $USAGE     => Help text

    my ($arg_href) = @_;

    ## Default(s)
    my $exit_code;

    ## Flatten argument(s)
    my $USAGE;

    my $tmpl = {
        exit_code => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$exit_code,
            strict_type => 1,
        },
        USAGE => {
            defined     => 1,
            required    => 1,
            store       => \$USAGE,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {*STDOUT} $USAGE;
    exit $exit_code;
}

sub create_temp_dir {

## Function : Create a temporary directory and returns the path to it
## Returns  : $temp_dir_path
## Arguments: $directory_base_name  => Base name of directroy
##          : $directory_path       => Where to create the temporary directory
##          : $FILEHANDLE           => Filehandle to write to
##          : $max_expression_value => Max integrer to add to directory_name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory_base_name;
    my $directory_path;
    my $FILEHANDLE;
    my $max_expression_value;

    my $tmpl = {
        directory_base_name => {
            default     => q{temp_dir},
            defined     => 1,
            store       => \$directory_base_name,
            strict_type => 1,
        },
        directory_path => {
            default     => cwd(),
            defined     => 1,
            store       => \$directory_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        max_expression_value => {
            default     => 1000,
            defined     => 1,
            store       => \$max_expression_value,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_mkdir };

    my $temp_dir_path = catdir( $directory_path,
        $directory_base_name . $UNDERSCORE . int rand $max_expression_value );

    if ($FILEHANDLE) {

        gnu_mkdir(
            {
                FILEHANDLE       => $FILEHANDLE,
                indirectory_path => $temp_dir_path,
                parents          => 1,
            }
        );
    }

    return $temp_dir_path;
}

sub print_parameter_defaults {

## Function : Print all parameters and their default values
## Returns  :
## Arguments: $colored                 => Colorize output
##          : $index                   => Display array indices
##          : $output                  => Output stream
##          : $parameter_href          => Holds all parameters {REF}
##          : $print_parameter_default => Print parameter default
##          : $scalar_quotes           => Quote symbols to enclose scalar values

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    ## Default
    my $colored;
    my $index;
    my $output;
    my $print_parameter_default;
    my $scalar_quotes;

    my $tmpl = {
        colored => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$colored,
            strict_type => 1,
        },
        index => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$index,
            strict_type => 1,
        },
        output => {
            default     => q{stderr},
            store       => \$output,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        print_parameter_default => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$print_parameter_default,
            strict_type => 1,
        },
        scalar_quotes => {
            allow       => [ undef, q{"} ],
            default     => q{"},
            store       => \$scalar_quotes,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Data::Printer;

    if ($print_parameter_default) {

        ## Print parameter hash an exit
        say {*STDERR} q{Default values from config file:};
        p(
            %{$parameter_href},
            colored       => $colored,
            index         => $index,
            scalar_quotes => $scalar_quotes,
            output        => $output,
        );

        exit 0;
    }
    return;
}

1;
