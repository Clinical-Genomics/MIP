package MIP::File::Path;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $DOT $FORWARD_SLASH $LOG_NAME $SINGLE_QUOTE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_allowed_temp_directory
      check_filesystem_objects_existance
      check_filesystem_objects_and_index_existance
      check_gzipped
      get_absolute_path
    };
}

sub check_allowed_temp_directory {

## Function : Check that the temp directory value is allowed
## Returns  :
## Arguments: $not_allowed_paths_ref => Not allowed paths for temp dir
##          : $temp_directory        => Temp directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $not_allowed_paths_ref;
    my $temp_directory;

    my $tmpl = {
        not_allowed_paths_ref => {
            default     => [],
            defined     => 1,
            store       => \$not_allowed_paths_ref,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %is_not_allowed = (
        $FORWARD_SLASH . q{scratch}                  => undef,
        $FORWARD_SLASH . q{scratch} . $FORWARD_SLASH => undef,
    );

    ## Add more than already defined paths
    map { $is_not_allowed{$_} = undef; } @{$not_allowed_paths_ref};

    # Return if value is allowed
    return 1 if ( not exists $is_not_allowed{$temp_directory} );

    $log->fatal( qq{$SINGLE_QUOTE--temp_directory }
          . $temp_directory
          . qq{$SINGLE_QUOTE is not allowed because MIP will remove the temp directory after processing.}
    );
    exit 1;
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
        ## Directory was found
        return 1 if ( -d $object_name );

        $error_msg =
          q{Could not find intended } . $parameter_name . q{ directory: } . $object_name;
        return ( 0, $error_msg );
    }
    ## Then object type must be file

    ## Check existence of supplied file
    ## File was found
    return 1 if ( -f $object_name );

    $error_msg =
      q{Could not find intended } . $parameter_name . q{ file: } . $object_name;
    return 0, $error_msg;
}

sub check_filesystem_objects_and_index_existance {

## Function : Checks if a file or directory file exists as well as index file. Croak if object or index file does not exist.
## Returns  :
## Arguments: $index_suffix   => Index file ending
##          : $is_build_file  => File object can be built
##          : $object_name    => Object to check for existance
##          : $object_type    => Type of item to check
##          : $parameter_name => MIP parameter name {REF}
##          : $path           => Path to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $object_name;
    my $object_type;
    my $parameter_name;
    my $path;

    ## Default
    my $index_suffix;
    my $is_build_file;

    my $tmpl = {
        index_suffix => {
            allow       => [qw{ gz .gz }],
            default     => q{.gz},
            store       => \$index_suffix,
            strict_type => 1,
        },
        is_build_file => {
            store       => \$is_build_file,
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
    return if ( defined $is_build_file );

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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

sub check_gzipped {

## Function : Check if a file is gzipped.
## Returns  : "0 (=uncompressed)| 1 (=compressed)"
## Arguments: $file_name => File name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## File gzipped
    return 1 if ( $file_name =~ / [.]gz$ /xms );

    ## File not gzipped
    return 0;
}

sub get_absolute_path {

## Function : Get absolute path for supplied path or croaks and exists if path does not exists
## Returns  : $path (absolute path)
## Arguments: $parameter_name => Parameter to be evaluated
##          : $path           => Supplied path to be updated/evaluated

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $parameter_name;
    my $path;

    my $tmpl = {
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Reformat to absolute path
    my $absolute_path = abs_path($path);

    return $absolute_path if ( defined $absolute_path );

    croak(  q{Could not find absolute path for }
          . $parameter_name . q{: }
          . $path
          . q{. Please check the supplied path!} );
}

1;
