#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Install::Post_installation} => [qw{ update_config }],
        q{MIP::Test::Fixtures}                      => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Install::Post_installation qw{ update_config };

diag(   q{Test update_config from Post_installation.pm v}
      . $MIP::Recipes::Install::Post_installation::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 0, } );

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given a config file path when it does not exists
my $env_name      = q{test};
my $pipeline      = q{rd_dna};
my $update_config = q{config_does_not_exists};
my $write_config  = 1;

trap {
    update_config(
        {
            env_name      => $env_name,
            filehandle    => $filehandle,
            pipeline      => $pipeline,
            update_config => $update_config,
            write_config  => $write_config,
        }
    )
};

## Then throw warning
like(
    $trap->stderr,
    qr/MIP\s+will\s+not\s+attempt\s+to\s+update/xms,
    q{Throw not update message}
);

## Close the filehandle
close $filehandle;

done_testing();
