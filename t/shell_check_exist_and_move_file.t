#!/usr/bin/env perl

use 5.026;
use Carp;
use Cwd;
use charnames qw{ :full :short };
use Digest::MD5;
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname  };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp qw{ tempfile };
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = 1.0.0;

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Language::Shell});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Language::Shell qw{ check_exist_and_move_file };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test check_exist_and_move_file from Shell.pm v}
      . $MIP::Language::Shell::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Set and create test directories
my $test_dir = catdir( $Bin, q{data} );
my $temp_dir = File::Temp->newdir();

## Create temporary filehandles
my $COMMAND_FILEHANDLE = File::Temp->new(
    TEMPLATE => q{James_XXXX},
    DIR      => $test_dir,
    SUFFIX   => q{.sh},
);
my $TEST_FILEHANDLE = File::Temp->new(
    TEMPLATE => q{Sir_Toby_XXXX},
    DIR      => $test_dir,
    SUFFIX   => q{.txt},
);
my $TEMP_FILEHANDLE_1 = File::Temp->new(
    TEMPLATE => q{Amiral_von_Scneider_XXXX},
    DIR      => $temp_dir,
    SUFFIX   => q{.txt},
);
my $TEMP_FILEHANDLE_2 = File::Temp->new(
    TEMPLATE => q{Mr_Pommeroy_XXXX},
    DIR      => $temp_dir,
    SUFFIX   => q{.txt},
);

## Write something to temporary files
say {$TEMP_FILEHANDLE_1} q{Cheerio Miss Sophie!};
say {$TEMP_FILEHANDLE_2} q{The same procedure as every year, James!};

## Get filenames
my $command_file_path = $COMMAND_FILEHANDLE->filename;
my $temp_file_path_1  = $TEMP_FILEHANDLE_1->filename;
my $temp_file_path_2  = $TEMP_FILEHANDLE_2->filename;
my $test_file_path    = $TEST_FILEHANDLE->filename;

## Create arrays to loop over
my @temp_paths       = ( $temp_file_path_1,  $temp_file_path_2 );
my @temp_filehandles = ( $TEMP_FILEHANDLE_1, $TEMP_FILEHANDLE_2 );

TEMP_PATH:
foreach my $temp_path (@temp_paths) {
    check_exist_and_move_file(
        {
            FILEHANDLE          => $COMMAND_FILEHANDLE,
            intended_file_path  => $test_file_path,
            temporary_file_path => $temp_path,
        }
    );
}

## Make files readable
my @filehandles =
  ( $COMMAND_FILEHANDLE, $TEMP_FILEHANDLE_1, $TEMP_FILEHANDLE_2, $TEST_FILEHANDLE );

FILEHANDLE:
foreach my $FILEHANDLE (@filehandles) {
    seek $FILEHANDLE, 0, 0
      or croak q{Seek on $COMMAND_FILEHANDLE failed:} . $ERRNO . $NEWLINE;
}

## Create MD5 object and md5 sum variables
my $md5 = Digest::MD5->new;
my $pre_test_file_md5;
my $post_test_file_md5;
my $temp_file_md5;

my $command_counter = 0;

COMMAND:
while ( my $command = <$COMMAND_FILEHANDLE> ) {
    next COMMAND if ( $command eq $NEWLINE );
    chomp $command;

    diag( q{Testing condition number} . $SPACE, $command_counter + 1 );

    ## Check that starting conditions are correct
    if ( $command_counter == 0 ) {
        ok( _file_is_zero( { file_path => $test_file_path } ) == 1,
            q{Test file starts empty} );
    }
    if ( $command_counter == 1 ) {
        ok( _file_has_size( { file_path => $test_file_path } ) == 1,
            q{Test file starts has size} );
    }
    ok( _file_has_size( { file_path => $temp_paths[$command_counter] } ) == 1,
        q{Temp file has size} );
    ## Capture md5 sum of files before command
    $pre_test_file_md5 = $md5->addfile($TEST_FILEHANDLE)->hexdigest;
    $temp_file_md5     = $md5->addfile( $temp_filehandles[$command_counter] )->hexdigest;
    ## Test that files start out as different
    isnt( $pre_test_file_md5, $temp_file_md5, q{Initial files have different md5sums} );

    ## Run command
    system $command;

    ## Test if the temporary file has been (re)moved
    ok( _file_exist( { file_path => $temp_paths[$command_counter] } ) == 0,
        q{Temp file has been removed} );

    ## Get new md5 sum on $test_file
    $post_test_file_md5 = $md5->addfile($TEST_FILEHANDLE)->hexdigest;
    if ( $command_counter == 0 ) {
        ok( _file_has_size( { file_path => $test_file_path } ) == 1,
            q{Test file has size after mv command} );
        is( $pre_test_file_md5, $post_test_file_md5,
            q{Test file has been replaced by temp file} );
    }
    if ( $command_counter == 1 ) {
        isnt( $post_test_file_md5, $temp_file_md5,
            q{Test file has not been replaced by temp file} );
    }

    $command_counter++;
}

close $COMMAND_FILEHANDLE;
close $TEMP_FILEHANDLE_1;
close $TEMP_FILEHANDLE_2;
close $TEST_FILEHANDLE;

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}

sub _file_is_zero {

## Function : Check if file has zero size
## Returns  : Returns 1 if true, otherwise 0
## Arguments: $file_path => Path to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;

    my $tmpl = {
        file_path => {
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Return 1 if file has size
    return 1 if ( -z $file_path );
    ## Otherwise return 0
    return 0;
}

sub _file_has_size {

## Function : Check if file has size
## Returns  : Returns 1 if true, otherwise 0
## Arguments: $file_path => Path to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;

    my $tmpl = {
        file_path => {
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Return 1 if file has size
    return 1 if ( -s $file_path );
    ## Otherwise return 0
    return 0;
}

sub _file_exist {

## Function : Check if file exists
## Returns  : Returns 1 if true, otherwise 0
## Arguments: $file_path => Path to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;

    my $tmpl = {
        file_path => {
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Return 1 if file exists
    return 1 if ( -e $file_path );
    ## Otherwise return 0
    return 0;
}
