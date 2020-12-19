package MIP::File::Path;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $DOT $EMPTY_STR $FORWARD_SLASH $LOG_NAME $SINGLE_QUOTE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_allowed_temp_directory
      check_filesystem_objects_existance
      check_filesystem_objects_and_index_existance
      get_absolute_path
      get_file_names
      get_file_line_by_line
      remove_file_path_suffix
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

    use MIP::Validate::Data qw{ %CONSTRAINT };

    my %exists_constraint_map = (
        directory => q{dir_exists},
        file      => q{plain_file_exists},
    );

    my $constraint = $exists_constraint_map{$object_type};
    return 1 if ( $CONSTRAINT{$constraint}->($object_name) );

    my $error_msg = qq{Could not find intended $parameter_name $object_type: $object_name};
    return ( 0, $error_msg );
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

sub get_file_names {

## Function : Get the file(s) from the filesystem
## Returns  : @file_names
## Arguments: $file_directory   => File directory
##          : $rule_name        => Rule name string
##          : $rule_skip_subdir => Rule skip sub directories

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_directory;
    my $rule_name;
    my $rule_skip_subdir;

    my $tmpl = {
        file_directory => {
            defined     => 1,
            required    => 1,
            store       => \$file_directory,
            strict_type => 1,
        },
        rule_name => {
            store       => \$rule_name,
            strict_type => 1,
        },
        rule_skip_subdir => {
            store       => \$rule_skip_subdir,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Path::Iterator::Rule;

    my @file_names;

    ## Get all files in supplied indirectories
    my $rule = Path::Iterator::Rule->new;

    ### Set rules
    ## Ignore if sub directory
    if ($rule_skip_subdir) {

        $rule->skip_subdirs($rule_skip_subdir);
    }

    ## Look for particular file name
    if ($rule_name) {

        $rule->name($rule_name);
    }

    # Initilize iterator
    my $iter = $rule->iter($file_directory);

  DIRECTORY:
    while ( my $file_path = $iter->() ) {

        my $file_name = splitpath($file_path);

        push @file_names, $file_name;
    }
    return @file_names;
}

sub get_file_line_by_line {

## Function  : Read file line by line and return array where each element is a line
## Returns   : \@lines
## Arguments : $chomp => Remove any end-of-line character sequences
##           : $path  => File path to read

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $path;

    ## Default(s)
    my $chomp;

    my $tmpl = {
        chomp => {
            default     => 0,
            store       => \$chomp,
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

    use Path::Tiny qw{ path };

    my @lines = path($path)->lines_utf8( { chomp => $chomp, } );
    return \@lines;
}

sub remove_file_path_suffix {

## Function : Parse file suffixes in file path. Removes suffix if matching else return undef
## Returns  : undef | $file_path_no_suffix
## Arguments: $file_path         => File path
##          : $file_suffixes_ref => File suffix to be removed

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $file_suffixes_ref;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        file_suffixes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$file_suffixes_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my ( $file_name_nosuffix, $dir, $suffix ) =
      fileparse( $file_path, @{$file_suffixes_ref} );

    $dir = $dir eq q{./} ? $EMPTY_STR : $dir;

    return $dir . $file_name_nosuffix if ( $file_name_nosuffix and $suffix );

    return;
}

1;
