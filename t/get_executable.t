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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $EMPTY_STR $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Environment::Executable} => [qw{ get_executable }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Executable qw{ get_executable };

diag(   q{Test get_executable from Executable.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an executable name
my %got_executable = get_executable( { executable_name => q{vep}, } );

my $expected_executable_href =
  { version_regexp =>
q?'my ($version) = /ensembl-vep\s+:\s(\d+)/xms; if($version) {print $version;last;}'?,
  };

## Then return hash ref for only executable
is_deeply( \%got_executable, $expected_executable_href, q{Got executable specific hash} );

## Given no executable name
my %got_all_executable = get_executable( {} );

my %expected_all_executable = (
    mip => {
        version_cmd => q{version},
        version_regexp =>
q?'my ($version) = /mip\s+version\s+(\S+)/xms; if($version) {print $version;last;}'?,
    },
);

## Then return hash ref for entire executable
is_deeply(
    $got_all_executable{mip},
    $expected_all_executable{mip},
    q{Got entire executable hash}
);

done_testing();
