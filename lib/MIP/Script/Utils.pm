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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      create_temp_dir
      help
      nest_hash
      print_install_defaults
      print_parameter_defaults
      update_program_versions
      write_script_version };
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
## Arguments: $parameter_href => Holds all parameters {REF}
##          : $colored        => Colorize output
##          : $index          => Display array indices
##          : $scalar_quotes  => Quote symbols to enclose scalar values

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    ## Default
    my $colored;
    my $index;
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
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

    ## Print parameter hash an exit
    say {*STDERR} q{Default values from config file:};
    p(
        %{$parameter_href},
        colored       => $colored,
        index         => $index,
        scalar_quotes => $scalar_quotes,
    );

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

    ## Nest the shell parameters
    my @installations = @{ $cmd_href->{installations} };
    my @colon_keys = grep { /:/xms } keys %{$cmd_href};
  PARAMETER:
    foreach my $parameter (@colon_keys) {

        my $final_value = $cmd_href->{$parameter};
        foreach my $installation (@installations) {
            _recursive_nesting(
                {
                    array_to_shift_ref => [ ( split /:/xms, $parameter ) ],
                    final_value => $final_value,
                    hash_to_populate_href => $cmd_href->{$installation},
                }
            );
        }
        delete $cmd_href->{$parameter};
    }
    return;
}

sub update_program_versions {
## Function : Set program versions in parameter hash
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}

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

    my @programs = keys %{ $parameter_href->{program_versions} };

    ## Check if any versions needs to be updated
    return if scalar @programs == 0;

  INSTALLATION:
    foreach my $installation ( @{ $parameter_href->{installations} } ) {

      PROGRAM:
        foreach my $program (@programs) {

          INSTALL_MODE:
            foreach my $install_mode ( keys %{ $parameter_href->{$installation} } ) {
                if ( $install_mode eq q{shell} ) {
                    if ( $parameter_href->{$installation}{$install_mode}{$program} ) {
                        $parameter_href->{$installation}{$install_mode}{$program}{version}
                          = $parameter_href->{program_versions}{$program};
                    }
                }
                else {
                    if ( $parameter_href->{$installation}{$install_mode}{$program} ) {
                        $parameter_href->{$installation}{$install_mode}{$program} =
                          $parameter_href->{program_versions}{$program};
                    }
                }
            }
        }
    }
    delete $parameter_href->{program_versions};
    return;
}

sub write_script_version {

## Function : Writes the vesrion of the script and exists
## Returns  :
## Arguments: $write_version => Write version
##          : $version       => Version of script to print

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $version;

    ## Default(s)
    my $write_version;

    my $tmpl = {
        write_version => {
            allow       => [ 0, 1, undef ],
            default     => 0,
            store       => \$write_version,
            strict_type => 1,
        },
        version => {
            defined     => 1,
            required    => 1,
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return 1 if ( not $write_version );

    say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $version, $NEWLINE;
    exit;
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
    my $value = shift @{$array_to_shift_ref};

    if ( not $hash_to_populate_href->{$value} ) {
        return;
    }

    ## If the array is empty, give the last hash key the final value and return
    if ( scalar @{$array_to_shift_ref} == 0 ) {
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
    return;
}

1;
