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
use MIP::Test::Fixtures qw{ test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Parameter}      => [qw{ set_cache_program_executables }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ set_cache_program_executables };

diag(   q{Test set_cache_program_executables from Parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given parameters with program executables and a cache recipe
my %parameter = test_mip_hashes( { mip_hash_name => q{define_parameter}, } );
push @{ $parameter{cache}{recipe} }, q{bwa_mem};

set_cache_program_executables( { parameter_href => \%parameter, } );

my %expected_cache = (
    cache => {
        program_executables => [qw{ bwa samtools sambamba }],
        recipe              => [qw{ bwa_mem }],
    },
);

## Then program executables should be set in cache
is_deeply( $parameter{cache}, $expected_cache{cache},
    q{Set program executables in cache} );

done_testing();
