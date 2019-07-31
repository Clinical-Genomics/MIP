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
use Modern::Perl qw{ 2017 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Parse::File}    => [qw{ parse_sing_bind_paths }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::File qw{ parse_sing_bind_paths };

diag(   q{Test parse_sing_bind_paths from Parse::File.pm v}
      . $MIP::Parse::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given paths with some ovelapping directories and duplicates
my @test_paths = (
    catdir(qw{ dir_a dir_c }),       catdir(qw{ dir_a dir_c }),
    catdir(qw{ dir_a dir_c dir_d }), catdir(qw{ dir_a }),
    catdir(qw{ dir_b dir_c }),       catdir(qw{ dir_b dir_c dir_d }),
    catdir(qw{ dir_b dir_e dir_d }), catdir(qw{ dir_d dir_a }),
);

my @return_paths = parse_sing_bind_paths( { dir_paths_ref => \@test_paths, } );

## Then remove duplicates and paths that with common roots
my @expected = (
    catdir(qw{ dir_a }),       catdir(qw{ dir_b dir_c }),
    catdir(qw{ dir_d dir_a }), catdir(qw{ dir_b dir_e dir_d }),
);

is_deeply( \@return_paths, \@expected, q{Return paths} );

done_testing();
