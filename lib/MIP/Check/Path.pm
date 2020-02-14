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
use List::MoreUtils qw{ any };
use Readonly;

## MIPs lib
use MIP::Constants qw{ $AMPERSAND $CLOSE_BRACKET $NEWLINE $OPEN_BRACKET $SPACE $TAB };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.12;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_file_version_exist
      check_future_filesystem_for_directory
      check_gatk_sample_map_paths
      check_target_bed_file_suffix
      check_vcfanno_toml
    };
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

sub check_future_filesystem_for_directory {

## Function : Build bash script to check if a directory exists and otherwise create it
## Returns  :
## Arguments: $directory_path => Path to check / create

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory_path;

    my $tmpl = {
        directory_path => {
            defined     => 1,
            required    => 1,
            store       => \$directory_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{gnu_mkdir};

    my $dir_check_command =
        $OPEN_BRACKET
      . $SPACE . q{! -d}
      . $SPACE
      . $directory_path
      . $SPACE
      . $CLOSE_BRACKET
      . $SPACE
      . $AMPERSAND
      . $AMPERSAND
      . $SPACE;
    my @mkdir_commands = gnu_mkdir(
        {
            indirectory_path => $directory_path,
            parents          => 1,
        }
    );

    $dir_check_command .= join $SPACE, @mkdir_commands;

    return $dir_check_command;
}

sub check_gatk_sample_map_paths {

## Function : Check that the supplied gatk sample map file paths exists
## Returns  :
## Arguments: $log             => Log object
##          : $sample_map_path => Sample map path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $sample_map_path;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        sample_map_path => {
            defined     => 1,
            required    => 1,
            store       => \$sample_map_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Constants
    Readonly my $RECORD_SEPARATOR => qq{\t};
    Readonly my $FIELD_COUNTER    => 2;

    ## Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    open $filehandle, q{<}, $sample_map_path
      or $log->logdie( q{Cannot open '} . $sample_map_path . q{': } . $OS_ERROR );

  LINE:
    while (<$filehandle>) {

        ## Remove newline
        chomp;

        ## Unpack line
        my $line = $_;

        ## Get sample and file path (and check proper format)
        my ( $sample, $file_path, $unexpected_data ) =
          split $RECORD_SEPARATOR, $line, $FIELD_COUNTER + 1;

        ## Make sure that we get what we expect
        if ( defined $unexpected_data ) {

            carp q{Unexpected trailing garbage at end of line '} . $line . q{':},
              $NEWLINE . $TAB . $unexpected_data . $NEWLINE;
        }
        ## Path exists
        next LINE if ( -e $file_path );

        $log->fatal( q{The supplied file path: }
              . $file_path
              . q{ from sample map file: }
              . $sample_map_path
              . q{ does not exist} );
        exit 1;
    }
    close $filehandle;
    return 1;
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

## Function : Check that the supplied vcfanno toml config has mandatory keys and file exists for annotation array
## Returns  :
## Arguments: $log               => Log object
##          : $parameter_name    => Name of parameter
##          : $vcfanno_file_toml => Toml config file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $parameter_name;
    my $vcfanno_file_toml;

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
        vcfanno_file_toml => {
            defined     => 1,
            required    => 1,
            store       => \$vcfanno_file_toml,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_and_index_existance };
    use MIP::Io::Read qw{ read_from_file };

    my %vcfanno_config = read_from_file(
        {
            path   => $vcfanno_file_toml,
            format => q{toml},
        }
    );
    my @vcfanno_features = qw{ file fields ops };
    my $err_msg = q{ is not defined or empty vcfanno toml features. Please check file: }
      . $vcfanno_file_toml;
  ANNOTATION:

    foreach my $annotation_href ( @{ $vcfanno_config{annotation} } ) {

      FEATURE:
        foreach my $feature (@vcfanno_features) {

            ## Check mandatory feature keys for vcfanno
            if ( not defined $annotation_href->{$feature} ) {

                $log->fatal( q{Feature: } . $feature . $err_msg );
                exit 1;
            }
        }

        ## Check path object exists
        check_filesystem_objects_and_index_existance(
            {
                object_name    => q{file},
                object_type    => q{file},
                parameter_name => $parameter_name,
                path           => $annotation_href->{file},
            }
        );
    }
    return 1;
}

1;
