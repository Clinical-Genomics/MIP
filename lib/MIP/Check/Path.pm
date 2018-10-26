package MIP::Check::Path;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use autodie;
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_command_in_path
      check_filesystem_objects_existance
      check_filesystem_objects_and_index_existance
      check_file_version_exist
      check_parameter_files
      check_target_bed_file_suffix
      check_vcfanno_toml
    };
}

## Constants
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};

sub check_command_in_path {

## Function : Checking commands in your path and executable
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
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

    use MIP::Check::Unix qw{check_binary_in_path};

    # Track program paths that have already been checked
    my %seen;

    ## Checking binaries for parameter program_source_environment_command
    %seen = _check_program_source_env_cmd_binary(
        {
            active_parameter_href => $active_parameter_href,
            log                   => $log,
        }
    );

    ## Checking program_name_path binaries
    _check_program_name_path_binaries(
        {
            active_parameter_href => $active_parameter_href,
            log                   => $log,
            parameter_href        => $parameter_href,
            seen_href             => \%seen,
        }
    );
    return;
}

sub check_filesystem_objects_existance {

## Function : Checks if a file or directory file exists
## Returns  : (0 | 1, $error_msg)
## Arguments: $object_name    => Object to check for existance
##          : $object_type    => Type of item to check
##          : $parameter_name => MIP parameter name {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $object_name;
    my $object_type;
    my $parameter_name;

    my $tmpl = {
        object_name => {
            defined     => 1,
            required    => 1,
            store       => \$object_name,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        object_type => {
            allow       => [qw{ directory file }],
            defined     => 1,
            required    => 1,
            store       => \$object_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## For potential error messages
    my $error_msg;

    ## Check existance of directory
    if ( $object_type eq q{directory} ) {

        ## Check existence of supplied directory
        if ( not -d $object_name ) {

            $error_msg =
                q{Could not find intended }
              . $parameter_name
              . q{ directory: }
              . $object_name;
            return ( 0, $error_msg );
        }
        ## File was found
        return 1;
    }
    ## Then object type must be file

    ## Check existence of supplied file
    if ( not -f $object_name ) {

        $error_msg =
          q{Could not find intended } . $parameter_name . q{ file: } . $object_name;

        return 0, $error_msg;
    }
    ## File was found
    return 1;
}

sub check_filesystem_objects_and_index_existance {

## Function : Checks if a file or directory file exists as well as index file. Croak if object or index file does not exist.
## Returns  :
## Arguments: $log            => Log object to write to
##          : $index_suffix   => Index file ending
##          : $object_name    => Object to check for existance
##          : $object_type    => Type of item to check
##          : $parameter_href => Parameters hash
##          : $parameter_name => MIP parameter name {REF}
##          : $path           => Path to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $object_name;
    my $object_type;
    my $parameter_href;
    my $parameter_name;
    my $path;

    ## Default
    my $index_suffix;

    my $tmpl = {
        log          => { required => 1, store => \$log, },
        index_suffix => {
            allow       => [qw{ .gz }],
            default     => q{gz},
            store       => \$index_suffix,
            strict_type => 1,
        },
        object_name => {
            defined     => 1,
            required    => 1,
            store       => \$object_name,
            strict_type => 1,
        },
        object_type => {
            allow       => [qw{ directory file }],
            defined     => 1,
            required    => 1,
            store       => \$object_type,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        path => {
            defined     => 1,
            required    => 1,
            store       => \$path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Special case for file with "build_file" in config
    ## These are handled downstream
    return if ( defined $parameter_href->{$parameter_name}{build_file} );

    my ( $exist, $error_msg ) = check_filesystem_objects_existance(
        {
            object_name    => $path,
            object_type    => $object_type,
            parameter_name => $parameter_name,
        }
    );
    if ( not $exist ) {
        $log->fatal($error_msg);
        exit 1;
    }

    ## Check for tabix index as well
    if ( $path =~ m{ $index_suffix$ }xsm ) {

        my $path_index = $path . $DOT . q{tbi};

        my ( $index_exist, $index_error_msg ) = check_filesystem_objects_existance(
            {
                object_name    => $path_index,
                object_type    => $object_type,
                parameter_name => $path_index,
            }
        );
        if ( not $index_exist ) {
            $log->fatal($index_error_msg);
            exit 1;
        }
    }
    return 1;
}

sub check_file_version_exist {

## Function : Check if a file with with a filename consisting of $file_path_prefix.$file_counter.$file_path_suffix exist. If so bumps the version number and return new file path and version number.
## Returns  : $file_path, $file_name_version
## Arguments: $file_path_prefix => File path prefix
##          : $file_path_suffix => File path suffix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path_prefix;
    my $file_path_suffix;

    my $tmpl = {
        file_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$file_path_prefix,
            strict_type => 1,
        },
        file_path_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$file_path_suffix,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Nr of sbatch scripts with identical filenames i.e. version number
    my $file_name_version = 0;

    my $file_path = $file_path_prefix . $file_name_version . $file_path_suffix;

  FILE_PATHS:
    while ( -e $file_path ) {

        $file_name_version++;

        ## New file_path to test for existence
        $file_path = $file_path_prefix . $file_name_version . $file_path_suffix;
    }
    return ( $file_path, $file_name_version );
}

sub check_parameter_files {

## Function : Checks that files/directories files exists
## Returns  :
## Arguments: $active_parameter_href   => Holds all set parameter for analysis
##          : $associated_recipes_ref  => The parameters recipe(s) {REF}
##          : $family_id               => The family_id
##          : $log                     => Log object to write to
##          : $parameter_exists_check  => Check if intendend file exists in reference directory
##          : $parameter_href          => Holds all parameters
##          : $parameter_name          => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $associated_recipes_ref;
    my $log;
    my $parameter_exists_check;
    my $parameter_href;
    my $parameter_name;

    ## Default(s)
    my $family_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        associated_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$associated_recipes_ref,
            strict_type => 1,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        log                    => { required => 1, store => \$log, },
        parameter_exists_check => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_exists_check,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Path
      qw{ check_filesystem_objects_existance check_filesystem_objects_and_index_existance };

    my %only_wgs = ( gatk_genotypegvcfs_ref_gvcf => 1, );

    ## Unpack parameters
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};

    ## Do nothing since parameter is not required unless exome mode is enabled
    return
      if ( exists $only_wgs{$parameter_name}
        && $consensus_analysis_type =~ / wgs /xsm );

    ## Check all recipes that use parameter
  ASSOCIATED_RECIPE:
    foreach my $associated_recipe ( @{$associated_recipes_ref} ) {

        ## Active associated recipe
        my $associated_recipe_name = $active_parameter_href->{$associated_recipe};

        ## Active parameter
        my $active_parameter = $active_parameter_href->{$parameter_name};

        ## Only check active associated recipes parameters
        next ASSOCIATED_RECIPE if ( not $associated_recipe_name );

        ## Only check active parameters
        next ASSOCIATED_RECIPE if ( not defined $active_parameter );

        if ( ref $active_parameter eq q{ARRAY} ) {

            ## Get path for array elements
          PATH:
            foreach my $path ( @{ $active_parameter_href->{$parameter_name} } ) {

                check_filesystem_objects_and_index_existance(
                    {
                        log            => $log,
                        object_name    => $path,
                        object_type    => $parameter_exists_check,
                        parameter_href => $parameter_href,
                        parameter_name => $parameter_name,
                        path           => $path,
                    }
                );
            }
            return;
        }
        elsif ( ref $active_parameter eq q{HASH} ) {

            ## Get path for hash keys
          PATH:
            for my $path ( keys %{ $active_parameter_href->{$parameter_name} } ) {

                check_filesystem_objects_and_index_existance(
                    {
                        log            => $log,
                        object_name    => $path,
                        object_type    => $parameter_exists_check,
                        parameter_href => $parameter_href,
                        parameter_name => $parameter_name,
                        path           => $path,
                    }
                );
            }
            return;
        }

        ## File
        my $path = $active_parameter_href->{$parameter_name};

        check_filesystem_objects_and_index_existance(
            {
                log            => $log,
                object_name    => $path,
                object_type    => $parameter_exists_check,
                parameter_href => $parameter_href,
                parameter_name => $parameter_name,
                path           => $path,
            }
        );
        return;
    }
    return;
}

sub check_target_bed_file_suffix {

## Function : Check that supplied target file ends with ".bed" and otherwise exists.
## Returns  :
## Arguments: $log            => Log object
##          : $parameter_name => MIP parameter name
##          : $path           => Path to check for ".bed" file ending

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $parameter_name;
    my $path;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $path !~ m{[.]bed$}xsm ) {

        $log->fatal(
            q{Could not find intendended '.bed file ending' for target file: }
              . $path
              . q{ in parameter '-}
              . $parameter_name . q{'},
            $NEWLINE
        );
        exit 1;
    }
    return 1;
}

sub check_vcfanno_toml {

## Function : Check that the supplied vcfanno toml frequency file match record 'file=' within toml config file
## Returns  :
## Arguments: $log               => Log object
##          : $vcfanno_file_freq => Frequency file recorded inside toml file
##          : $vcfanno_file_toml => Toml config file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $vcfanno_file_freq;
    my $vcfanno_file_toml;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        vcfanno_file_freq => {
            defined     => 1,
            required    => 1,
            store       => \$vcfanno_file_freq,
            strict_type => 1,
        },
        vcfanno_file_toml => {
            defined     => 1,
            required    => 1,
            store       => \$vcfanno_file_toml,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    open $FILEHANDLE, q{<}, $vcfanno_file_toml
      or $log->logdie( q{Cannot open '} . $vcfanno_file_toml . q{': } . $OS_ERROR );

  LINE:
    while (<$FILEHANDLE>) {

        ## Remove newline
        chomp;

        my ($file_path_freq) = /^file="(\S+)"/sxm;

        next LINE if ( not $file_path_freq );

        if ( $file_path_freq ne $vcfanno_file_freq ) {

            $log->fatal( q{The supplied vcfanno_config_file: }
                  . $vcfanno_file_freq
                  . q{ does not match record 'file=}
                  . $file_path_freq
                  . q{' in the sv_vcfanno_config file: }
                  . $vcfanno_file_toml );
            exit 1;
        }
        else {
            close $FILEHANDLE;
            return 1;
        }
    }
    return;
}

sub _check_program_name_path_binaries {

## Function : Checking program_name_path binaries
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $parameter_href        => Parameter hash {REF}
##          : $seen_href             => Track program paths already checked {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $parameter_href;
    my $seen_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        seen_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$seen_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Unix qw{check_binary_in_path};

  PARAMETER:
    foreach my $parameter_name ( keys %{$active_parameter_href} ) {

        ## Only check path(s) for parameters with "type" key
        next PARAMETER
          if ( not exists $parameter_href->{$parameter_name}{type} );

        ## Only check path(s) for parameters with type value eq "recipe"
        next PARAMETER
          if ( not $parameter_href->{$parameter_name}{type} eq q{recipe} );
        ## Only check path(s) for active recipes
        next PARAMETER if ( not $active_parameter_href->{$parameter_name} );

        ## Alias
        my $program_name_paths_ref =
          \@{ $parameter_href->{$parameter_name}{program_name_path} };

      PROGRAM:
        foreach my $program ( @{$program_name_paths_ref} ) {

            ## Only check path once
            next PROGRAM if ( $seen_href->{$program} );

            $seen_href->{$program} = check_binary_in_path(
                {
                    active_parameter_href => $active_parameter_href,
                    binary                => $program,
                    log                   => $log,
                    program_name          => $parameter_name,
                }
            );
        }
    }
    return;
}

sub _check_program_source_env_cmd_binary {

## Function : Checking commands for program source environment command binary exists and is executable
## Returns  : %seen
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Unix qw{check_binary_in_path};

    # Track program paths that have already been checked
    my %seen;

    foreach my $program (
        keys %{ $active_parameter_href->{program_source_environment_command} } )
    {

        ## Program program source env cmd binaries
        $seen{$program} = check_binary_in_path(
            {
                active_parameter_href => $active_parameter_href,
                binary                => $program,
                log                   => $log,
                program_name          => $program,
            }
        );
    }
    return %seen;
}

1;
