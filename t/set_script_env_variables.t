#! /usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Script::Setup_script} => [qw{ set_script_env_variables }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Setup_script qw{ set_script_env_variables };

diag(   q{Test set_script_env_variables from Setup_script.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given a temp dir and xdg_runtime_dir
my $temp_directory_bash = q{a_temp_dir};
my $xdg_runtime_dir     = 1;

## When setting env in script
set_script_env_variables(
    {
        filehandle          => $filehandle,
        temp_directory_bash => $temp_directory_bash,
        xdg_runtime_dir     => $xdg_runtime_dir,
    }
);

## Close the filehandle
close $filehandle;

my $expected_xdg_runtime_dir = q{XDG_RUNTIME_DIR=} . $temp_directory_bash;

## Then the env variable should be written
my ($returned_base_command) = $file_content =~ /^($expected_xdg_runtime_dir)/ms;

ok( $returned_base_command, q{Wrote xdg_runtime_dir to file} );
done_testing();
