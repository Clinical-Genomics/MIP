#!/usr/bin/env perl

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
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Environment::Executable} =>
          [qw{ build_binary_version_cmd get_executable }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Executable qw{ build_binary_version_cmd get_executable };

diag(   q{Test build_binary_version_cmd from Executable.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a binary
my $binary     = q{genmod};
my %executable = get_executable( { executable_name => $binary, } );

my @version_cmds = build_binary_version_cmd(
    {
        binary_path    => $binary,
        version_cmd    => $executable{version_cmd},
        version_regexp => $executable{version_regexp},
    }
);
my $version_cmds = join $SPACE, @version_cmds;

## Then return should contain version cmd
like( $version_cmds, qr/genmod\s+--version/xms, q{Got version cmd} );

## Then return should contain version regexp
like( $version_cmds, qr/version;last;/xms, q{Got version regexp} );

done_testing();
