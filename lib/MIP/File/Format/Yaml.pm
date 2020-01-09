package MIP::File::Format::Yaml;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Readonly;
use YAML qw{ Dump LoadFile };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ load_yaml order_parameter_names write_yaml };
}

## Constants
Readonly my $COLON        => q{:};
Readonly my $DOT          => q{.};
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };
Readonly my $NEWLINE      => qq{\n};

sub load_yaml {

## Function : Loads a YAML file into an arbitrary hash and returns it.
## Returns  : %yaml
## Arguments: $yaml_file => The yaml file to load

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $yaml_file;

    my $tmpl = {
        yaml_file => {
            defined     => 1,
            required    => 1,
            store       => \$yaml_file,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %yaml;

    open my $YAML, q{<}, $yaml_file
      or croak q{cannot open} . $SPACE . $yaml_file . $COLON . $OS_ERROR,
      $NEWLINE;

    ## Load hashreference as hash
    %yaml = %{ LoadFile($yaml_file) };

    close $YAML;

    return %yaml;
}

sub write_yaml {

## Function : Writes a YAML hash to file
## Returns  :
## Arguments: $yaml_file_path => Yaml file to write to
##          : $yaml_href      => Hash to dump {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $yaml_file_path;
    my $yaml_href;

    my $tmpl = {
        yaml_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$yaml_file_path,
            strict_type => 1,
        },
        yaml_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$yaml_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Localize the current value
    local $YAML::QuoteNumericStrings = 1;

    open my $YAML, q{>}, $yaml_file_path
      or croak q{Cannot open}
      . $SPACE
      . $SINGLE_QUOTE
      . $DOUBLE_QUOTE
      . $yaml_file_path
      . $DOUBLE_QUOTE
      . $SINGLE_QUOTE
      . $COLON
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
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Filehandles
    my $filehandle = IO::Handle->new();

    open $filehandle, q{<}, $file_path
      or croak( q{Cannot open}
          . $DOT
          . $SINGLE_QUOTE
          . $file_path
          . $SINGLE_QUOTE
          . $COLON
          . $SPACE
          . $OS_ERROR
          . $NEWLINE );

    my @order_parameters = _parse_yaml_file( { filehandle => $filehandle, } );

    close $filehandle;
    return @order_parameters;
}

sub _parse_yaml_file {

## Returns  : @order_parameters
## Arguments: $filehandle  => Filehandle to read

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;

    my $tmpl = { filehandle => { store => \$filehandle, }, };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Hold the order of the parameters from yaml file
    my @order_parameters;

  LINE:
    while ( my $line = <$filehandle> ) {

        chomp $line;

        ## Next line if header
        next LINE if ( $INPUT_LINE_NUMBER == 1 && $line =~ /---/sxm );

        ## Next line if commment
        next LINE if ( $line =~ /^#/sm );

        ## First level key
        if ( $line =~ /^(\w+):/sxm ) {

            my $parameter_name = $1;

            ## Add to enable later evaluation of parameters in proper order & write to MIP log file
            push @order_parameters, $parameter_name;
            next LINE;
        }
    }
    return @order_parameters;
}

1;
