package MIP::Check::Path;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ check_filesystem_objects_existance check_filesystem_objects_and_index_existance check_file_version_exist check_dir_path_exist check_target_bed_file_exist check_parameter_files };
}

## Constants
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};

sub check_filesystem_objects_and_index_existance {

## Function : Checks if a file or directory file exists as well as index file. Croak if object or index file does not exist.
## Returns  :
## Arguments: $parameter_href => Parameters hash
##          : $object_name    => Object to check for existance
##          : $parameter_name => MIP parameter name {REF}
##          : $object_type    => Type of item to check
##          : $path           => Path to check
##          : $index_suffix   => Index file ending

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $object_name;
    my $parameter_name;
    my $object_type;
    my $path;
    my $index_suffix;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        object_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$object_name,
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
        object_type => {
            required    => 1,
            defined     => 1,
            allow       => [qw{ directory file }],
            strict_type => 1,
            store       => \$object_type,
        },
        path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$path,
        },
        index_suffix => {
            default     => q{gz},
            allow       => [qw{ .gz }],
            strict_type => 1,
            store       => \$index_suffix
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Special for file with "build_file" in config
    ## These are handled downstream
    return if ( defined $parameter_href->{$parameter_name}{build_file} );

    my ( $exist, $error_msg ) = check_filesystem_objects_existance(
        {
            object_name    => $path,
            parameter_name => $parameter_name,
            object_type    => $object_type,
        }
    );
    if ( not $exist ) {
        $log->fatal($error_msg);
        exit 1;
    }

    ## Check for tabix index as well
    if ( $path =~ m{ $index_suffix$ }xsm ) {

        my $path_index = $path . $DOT . q{tbi};

        my ( $index_exist, $index_error_msg ) =
          check_filesystem_objects_existance(
            {
                object_name    => $path_index,
                parameter_name => $path_index,
                object_type    => $object_type,
            }
          );
        if ( not $index_exist ) {
            $log->fatal($index_error_msg);
            exit 1;
        }
    }
    return;
}

sub check_filesystem_objects_existance {

## Function : Checks if a file or directory file exists
## Returns  : (0 | 1, $error_msg)
## Arguments: $object_name    => Object to check for existance
##          : $parameter_name => MIP parameter name {REF}
##          : $object_type    => Type of item to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $object_name;
    my $parameter_name;
    my $object_type;

    my $tmpl = {
        object_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$object_name,
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
        object_type => {
            required    => 1,
            defined     => 1,
            allow       => [qw{ directory file }],
            strict_type => 1,
            store       => \$object_type,
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
        return 1;
    }
    elsif ( $object_type eq q{file} ) {

        ## Check existence of supplied file
        if ( not -f $object_name ) {

            $error_msg =
                q{Could not find intended }
              . $parameter_name
              . q{ file: }
              . $object_name;

            return 0, $error_msg;
        }
        else {
            ## File was found
            return 1;
        }
    }
    return;
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
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path_prefix,
        },
        file_path_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path_suffix,
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

sub check_dir_path_exist {

## Function  : Checks if any of the supplied paths exists. Returns an array that holds the existing paths.
## Returns   : @existing_dir_paths
## Arguments : $dir_paths_ref => Ref to array of supplied paths

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dir_paths_ref;

    my $tmpl = {
        dir_paths_ref => {
            default     => [],
            required    => 1,
            strict_type => 1,
            store       => \$dir_paths_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @existing_dir_paths;

    ## Search for directory in supplied paths
    @existing_dir_paths = grep { -d } @{$dir_paths_ref};

    return @existing_dir_paths;
}

sub check_target_bed_file_exist {

## Function : Check that supplied target file ends with ".bed" and otherwise exists.
## Returns  :
## Arguments: $parameter_name => MIP parameter name
##            $path           => Path to check for ".bed" file ending

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_name;
    my $path;

    my $tmpl = {
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

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

sub check_parameter_files {

## Function : Checks that files/directories files exists
## Returns  :
## Arguments: $parameter_href          => Holds all parameters
##          : $active_parameter_href   => Holds all set parameter for analysis
##          : $associated_programs_ref => The parameters program(s) {REF}
##          : $family_id               => The family_id
##          : $parameter_name          => Parameter name
##          : $parameter_exists_check  => Check if intendent file exists in reference directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $associated_programs_ref;
    my $parameter_name;
    my $parameter_exists_check;

    ## Default(s)
    my $family_id;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        associated_programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$associated_programs_ref
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
        parameter_exists_check => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_exists_check
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Path
      qw{ check_filesystem_objects_existance check_filesystem_objects_and_index_existance };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my %only_wgs = ( gatk_genotypegvcfs_ref_gvcf => 1, );

    ## Unpack parameters
    my $consensus_analysis_type =
      $parameter_href->{dynamic_parameter}{consensus_analysis_type};

    ## Do nothing since parameter is not required unless exome mode is enabled
    return
      if ( exists $only_wgs{$parameter_name}
        && $consensus_analysis_type =~ / wgs /xsm );

    ## Check all programs that use parameter
  ASSOCIATED_PROGRAM:
    foreach my $associated_program ( @{$associated_programs_ref} ) {

        ## Active associated program
        my $associated_program_name =
          $active_parameter_href->{$associated_program};

        ## Active parameter
        my $active_parameter = $active_parameter_href->{$parameter_name};

        ## Only check active associated programs parameters
        next ASSOCIATED_PROGRAM if ( not $associated_program_name );

        ## Only check active parameters
        next ASSOCIATED_PROGRAM if ( not defined $active_parameter );

        if ( ref $active_parameter eq q{ARRAY} ) {

            ## Get path for array elements
          PATH:
            foreach my $path ( @{ $active_parameter_href->{$parameter_name} } )
            {

                check_filesystem_objects_and_index_existance(
                    {
                        parameter_href => $parameter_href,
                        object_name    => $path,
                        parameter_name => $parameter_name,
                        object_type    => $parameter_exists_check,
                        path           => $path,
                    }
                );
            }
            return;
        }
        elsif ( ref $active_parameter eq q{HASH} ) {

            ## Get path for hash keys
          PATH:
            for my $path ( keys %{ $active_parameter_href->{$parameter_name} } )
            {

                check_filesystem_objects_and_index_existance(
                    {
                        parameter_href => $parameter_href,
                        object_name    => $path,
                        parameter_name => $parameter_name,
                        object_type    => $parameter_exists_check,
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
                parameter_href => $parameter_href,
                object_name    => $path,
                parameter_name => $parameter_name,
                object_type    => $parameter_exists_check,
                path           => $path,
            }
        );
        return;
    }
    return;
}

1;
