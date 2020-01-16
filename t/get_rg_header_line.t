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
        q{MIP::Sample_info}    => [qw{ get_rg_header_line }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ get_rg_header_line };

diag(   q{Test get_rg_header_line from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given input parameters
my %sample_info = (
    sample => {
        ADM1059A1 => {
            file => {
                ADM1059A1_161011_TestFilev2_GAGATTCC_lane1 => {
                    read_direction_file => {
                        ADM1059A1_161011_TestFilev2_GAGATTCC_lane1_1 => {
                            flowcell       => q{TestFilev2},
                            lane           => q{1},
                            sample_barcode => q{GAGATTC},
                        },
                    },
                },
            },
        },
    },
);

my %active_parameter = ( platform => q{ILLUMINA}, );
my $infile_prefix    = q{ADM1059A1_161011_TestFilev2_GAGATTCC_lane1};

my $rg_header_line = get_rg_header_line(
    {
        infile_prefix    => $infile_prefix,
        platform         => $active_parameter{platform},
        sample_id        => q{ADM1059A1},
        sample_info_href => \%sample_info,
        separator        => q{\t},
    }
);

## Then return expected read group header line
my $expected_rg_header_line =
q{ID:ADM1059A1_161011_TestFilev2_GAGATTCC_lane1\tLB:ADM1059A1\tPL:ILLUMINA\tPU:TestFilev2.1.GAGATTC\tSM:ADM1059A1};
is( $rg_header_line, $expected_rg_header_line, q{Get read group header line} );

done_testing();
