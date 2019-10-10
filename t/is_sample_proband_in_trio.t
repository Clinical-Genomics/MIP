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
        q{MIP::File::Format::Pedigree} => [qw{ is_sample_proband_in_trio }],
        q{MIP::Test::Fixtures}         => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Pedigree qw{ is_sample_proband_in_trio };

diag(   q{Test is_sample_proband_in_trio from Pedigre.pm v}
      . $MIP::File::Format::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %sample_info = test_mip_hashes( { mip_hash_name => q{qc_sample_info}, } );

## Given a non trio
$sample_info{has_trio} = 0;
my $check = is_sample_proband_in_trio(
    {
        sample_id        => q{ADM1059A1},
        sample_info_href => \%sample_info,
    }
);
## Then return 0
is( $check, 0, q{Detect not trio} );

## Given a unaffected
$sample_info{has_trio} = 1;
$check = is_sample_proband_in_trio(
    {
        sample_id        => q{ADM1059A2},
        sample_info_href => \%sample_info,
    }
);
## Then return 0
is( $check, 0, q{Not affected} );

## Given a affected child in trio
$check = is_sample_proband_in_trio(
    {
        sample_id        => q{ADM1059A1},
        sample_info_href => \%sample_info,
    }
);
## Then return 1
ok( $check, q{Detect affected child in trio} );

## Given an affected that is not a child
$sample_info{sample}{ADM1059A12}{phenotype} = q{affected};
$check = is_sample_proband_in_trio(
    {
        sample_id        => q{ADM1059A2},
        sample_info_href => \%sample_info,
    }
);
## Then return 0
is( $check, 0, q{Not child} );

done_testing();
