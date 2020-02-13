#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ remove_tree };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $JAVA_MEMORY_ALLOCATION => 4;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Xargs} => [qw{ xargs_command }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Xargs qw{ xargs_command };

diag(   q{Test xargs_command from Xargs.pm v}
      . $MIP::Recipes::Analysis::Xargs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

# Create anonymous filehandle
my $filehandle      = IO::Handle->new();
my $xargsfilehandle = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given analysis parameters
my $recipe_name            = q{xargs};
my $xargs_file_counter     = 0;
my $xargs_file_path_prefix = q{a_prefix};

## Create file commands for xargs
( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
    {
        core_number          => 1,
        filehandle           => $filehandle,
        file_path            => q{a_file_path},
        first_command        => q{bcftools},
        java_jar             => catfile(q{picard.jar}),
        java_use_large_pages => 1,
        memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
        recipe_info_path     => q{a_file_path},
        temp_directory       => q{a_temp_dir},
        xargsfilehandle      => $xargsfilehandle,
        xargs_file_counter   => $xargs_file_counter,
    }
);

## Create file commands for xargs
( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
    {
        core_number          => 1,
        filehandle           => $filehandle,
        file_path            => q{a_file_path},
        first_command        => q{java},
        java_jar             => catfile(q{picard.jar}),
        java_use_large_pages => 1,
        memory_allocation    => q{Xmx} . $JAVA_MEMORY_ALLOCATION . q{g},
        recipe_info_path     => q{a_file_path},
        temp_directory       => q{a_temp_dir},
        xargsfilehandle      => $xargsfilehandle,
        xargs_file_counter   => $xargs_file_counter,
    }
);

## Close the filehandle
close $filehandle;

## Then return TRUE
ok( $xargs_file_counter, q{ Executed analysis recipe } . $recipe_name );

## Clean-up
remove_tree(q{a_file_path.0.xargs});
remove_tree(q{a_file_path.1.xargs});

done_testing();
