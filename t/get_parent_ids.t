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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Qccollect}      => [qw{ get_parent_ids }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ get_parent_ids };

diag(   q{Test get_parent_ids from Qccollect.pm v}
      . $MIP::Qccollect::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a case with sample id and a pedigree
my %case = (
    ADM1059A1 => undef,
    ADM1059A2 => undef,
    ADM1059A3 => undef,
);

my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => q{a_recipe},
    }
);

my ( $father_id, $mother_id ) = get_parent_ids(
    {
        case_href        => \%case,
        sample_info_href => \%sample_info,
    }
);

## Then return father and mother id
is( $father_id, q{ADM1059A2}, q{Got father id} );
is( $mother_id, q{ADM1059A3}, q{Got mother id} );

## Given no parents
$sample_info{sample}{ADM1059A1}{father} = 0;
$sample_info{sample}{ADM1059A1}{mother} = 0;

( $father_id, $mother_id ) = get_parent_ids(
    {
        case_href        => \%case,
        sample_info_href => \%sample_info,
    }
);

## Then return fake father and mother id
is( $father_id, q{YYY_father}, q{Got fake father id} );
is( $mother_id, q{XXX_mother}, q{Got fake mother id} );

done_testing();
