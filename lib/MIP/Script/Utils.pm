package MIP::Script::Utils;

use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ 
        create_temp_dir 
        help 
        nest_hash 
        print_install_defaults 
        set_default_array_parameters };
}

## Constants
Readonly my $COLON      => q{:};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq(\n);
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

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
            defined => 1, 
            required => 1, 
            store => \$USAGE, 
            strict_type => 1, 
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {*STDOUT} $USAGE;
    exit $exit_code;
}

sub set_default_array_parameters {

## Function : Set default for array parameters unless parameter already exists in parameter hash
## Returns  :
## Arguments: $array_parameter_href => Hold the array parameter defaults as {REF}
##          : $parameter_href       => Parameters hash {REF}


    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $array_parameter_href;

    my $tmpl = {
        array_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$array_parameter_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
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
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
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

    gnu_mkdir(
        {
            FILEHANDLE       => $FILEHANDLE,
            indirectory_path => $temp_dir_path,
            parents          => 1,
        }
    );

    return $temp_dir_path;
}

sub print_install_defaults {

## Function : Print all parameters and their default values
## Returns  :
## Arguments: $parameter_href => Holds all parameters {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set default for vep cache dir
    if ( $parameter_href->{shell}{vep} ) {
        $parameter_href->{shell}{vep}{vep_cache_dir} = catdir( qw{ PATH TO CONDA },
            q{ensembl-tools-release-} . $parameter_href->{shell}{vep}{version},
            q{cache} );
    }

    ## Looping over the parameter hash to extract keys and values
  KEY:
    foreach my $key ( keys %{$parameter_href} ) {
        ## If the first level value is not a hash or array ref
        if ( ref( $parameter_href->{$key} ) !~ / ARRAY | HASH /xms ) {
            print {*STDOUT} $key . $SPACE;
            ## Check if scalar exists and print
            if ( $parameter_href->{$key} ) {
                say {*STDOUT} $parameter_href->{$key};
            }
            ## Boolean value
            else {
                say {*STDOUT} q{0};
            }
        }
        ## If the first level value is a hash ref
        elsif ( ref( $parameter_href->{$key} ) =~ /HASH/xms ) {
            ## Loop over the next set of hash keys
          PROGRAM:
            foreach my $program ( keys %{ $parameter_href->{$key} } ) {
                ## If the value is a hash ref
                if ( ref( $parameter_href->{$key}{$program} ) =~ /HASH/xms ) {
                    ## Loop over the next set of hash keys
                  NESTED_PARAM:
                    foreach my $nested_param (
                        keys %{ $parameter_href->{$key}{$program} } )
                    {
                        ## Print the key
                        print {*STDOUT} $key
                          . $SPACE
                          . $program
                          . $SPACE
                          . $nested_param
                          . $COLON
                          . $SPACE;
                        ## If the value is an array ref
                        if (
                            ref(
                                $parameter_href->{$key}{$program}{$nested_param}
                            ) =~ /ARRAY/xms
                          )
                        {
                            ## Print array
                            say {*STDOUT} join $SPACE,
                              @{ $parameter_href->{$key}{$program}
                                  {$nested_param} };
                        }
                        else {
                            ## Otherwise print the hash value
                            say {*STDOUT}
                              $parameter_href->{$key}{$program}{$nested_param};
                        }
                    }
                }
                ## Print values
                else {
                    ## Don't print value if it is undef
                    if ( not $parameter_href->{$key}{$program} ) {
                        say {*STDOUT} $key . $SPACE . $program;
                    }
                    else {
                        ## Print hash value
                        say {*STDOUT} $key
                          . $SPACE
                          . $program
                          . $COLON
                          . $SPACE
                          . $parameter_href->{$key}{$program};
                    }
                }
            }
        }
        ## Check for ref to array and print
        elsif ( ref( $parameter_href->{$key} ) =~ /ARRAY/xms ) {
            say {*STDOUT} $key . $COLON . $SPACE . join $SPACE,
              @{ $parameter_href->{$key} };
        }
    }
    exit 0;
}

sub nest_hash {
## Function : If necessary, nests the command line hash to fit the structure used in mip_install. 
##          : Splits the key string on ":" and creates nested keys from the split.
## Returns  :
## Arguments: $cmd_ref => Arguments from command line {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cmd_href;

    my $tmpl = {
        cmd_href => {
            defined  => 1,
            required => 1,
            store    => \$cmd_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Fix the shell programs into the hash
    if ( $cmd_href->{shell_programs} ) {
      PROGRAM:
        foreach my $program ( keys %{ $cmd_href->{shell_programs} } ) {

            $cmd_href->{shell}{$program}{version} =
              $cmd_href->{shell_programs}{$program};
            delete $cmd_href->{shell_programs}{$program};
        }
    }

    ## Nest the shell parameters
    my @colon_keys = grep { /:/ } keys %{$cmd_href};
  PARAMETER:
    foreach my $parameter (@colon_keys) {

        my $final_value = $cmd_href->{$parameter};
        _recursive_nesting(
            {
                array_to_shift_ref => [ ( split /:/, $parameter ) ],
                final_value => $final_value,
                hash_to_populate_href => $cmd_href,
            }
        );
        delete $cmd_href->{$parameter};
    }
    return;
}

sub _recursive_nesting {

## Function  : Recursive sub to nest values into a hash from an array
## Returns   :
## Arguments : $array_to_shift_ref    => Array of keys
##           : $final_value           => Value to be stored
##           : $hash_to_populate_href => Shift array values to this hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $array_to_shift_ref;
    my $final_value;
    my $hash_to_populate_href;

    my $tmpl = {
        final_value => {
            required => 1,
            store    => \$final_value,
        },
        hash_to_populate_href => {
            defined  => 1,
            required => 1,
            store    => \$hash_to_populate_href,
        },
        array_to_shift_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$array_to_shift_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Assign and remove the first value from the array
    my $value = shift @{ $array_to_shift_ref };

    ## If the array is empty, give the last hash key the final value and return
    if ( scalar @{ $array_to_shift_ref } == 0 ) {
        return $hash_to_populate_href->{$value} = $final_value;
    }

    ## Call same subroutine but increment the hash_ref to include the value as a key
    _recursive_nesting(
        {
            hash_to_populate_href => \%{ $hash_to_populate_href->{$value} },
            final_value           => $final_value,
            array_to_shift_ref    => $array_to_shift_ref,
        }
    );
}

1;
