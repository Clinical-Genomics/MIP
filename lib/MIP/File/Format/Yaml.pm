package MIP::File::Format::Yaml;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;
use YAML;
use autodie;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ load_yaml write_yaml order_parameter_names };
}

## Constants
Readonly my $DOT          => q{.};
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };
Readonly my $NEWLINE      => qq{\n};

# Localize the current value
local $YAML::QuoteNumericStrings = $YAML::QuoteNumericStrings;

# Force numeric values to strings in YAML representation
$YAML::QuoteNumericStrings = 1;

sub load_yaml {

## Function : Loads a YAML file into an arbitrary hash and returns it.
## Returns  : %yaml
## Arguments: $yaml_file => The yaml file to load

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $yaml_file;

    my $tmpl = {
        yaml_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$yaml_file,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %yaml;

    open my $YAML, q{<}, $yaml_file
      or croak q{cannot open} . $SPACE . $yaml_file . q{:} . $OS_ERROR,
      $NEWLINE;

    ## Load hashreference as hash
    %yaml = %{ YAML::LoadFile($yaml_file) };

    close $YAML;

    return %yaml;
}

sub write_yaml {

## Function : Writes a YAML hash to file
## Returns  :
## Arguments: $yaml_href      => The hash to dump {REF}
##          : $yaml_file_path => The yaml file to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $yaml_href;
    my $yaml_file_path;

    my $tmpl = {
        yaml_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$yaml_href,
        },
        yaml_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$yaml_file_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    open my $YAML, q{>}, $yaml_file_path
      or croak q{Cannot open}
      . $SPACE
      . $SINGLE_QUOTE
      . $DOUBLE_QUOTE
      . $yaml_file_path
      . $DOUBLE_QUOTE
      . $SINGLE_QUOTE . q{:}
      . $OS_ERROR
      . $NEWLINE;

    say {$YAML} Dump($yaml_href);
    close $YAML;
    return;
}

sub order_parameter_names {

##Function : Adds the order of first level keys from yaml file to array
##Returns  : @order_parameters
##Arguments: $file_path => File path to yaml file

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;

    my $tmpl = {
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Filehandles
    my $FILEHANDLE = IO::Handle->new();

    open $FILEHANDLE, q{<}, $file_path
      or croak( q{Cannot open}
          . $DOT
          . $SINGLE_QUOTE
          . $file_path
          . $SINGLE_QUOTE . q{:}
          . $SPACE
          . $OS_ERROR
          . $NEWLINE );

    my @order_parameters = _parse_yaml_file( { filehandle => $FILEHANDLE, } );

    close $FILEHANDLE;
    return @order_parameters;
}

sub _parse_yaml_file {

## Returns  : @order_parameters
## Arguments: $FILEHANDLE  => Filehandle to read

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = { filehandle => { store => \$FILEHANDLE }, };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Hold the order of the parameters from yaml file
    my @order_parameters;

  LINE:
    while (<$FILEHANDLE>) {

        chomp;

        ## Next line if header
        next LINE if ( $INPUT_LINE_NUMBER == 1 && m/---/sxm );

        ## Next line if commment
        next LINE if (/^#/sm);

        ## First level key
        if (m/^(\w+):/sxm) {

            my $parameter_name = $1;

            ## Add to enable later evaluation of parameters in proper order & write to MIP log file
            push @order_parameters, $parameter_name;
            next LINE;
        }
    }
    return @order_parameters;
}

1;
