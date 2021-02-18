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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ get_pedigree_sample_ids }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ get_pedigree_sample_ids };

diag(   q{Test get_pedigree_sample_ids from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given sample ids in sample_info
my %sample_info = test_mip_hashes( { mip_hash_name => q{qc_sample_info} } );

## When getting sample ids
my @sample_ids = get_pedigree_sample_ids( { sample_info_href => \%sample_info } );

## Sorting to make test predictable
@sample_ids = sort @sample_ids;

my @expected_sample_ids = sort qw{ ADM1059A1 ADM1059A2 ADM1059A3 };

## Then return samples ids
is_deeply( \@sample_ids, \@expected_sample_ids, q{Got sample ids} );

done_testing();
