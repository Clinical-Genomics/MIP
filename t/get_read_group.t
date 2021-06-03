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


## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ get_read_group }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ get_read_group };

diag(   q{Test get_read_group from Sample_info.pm}
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

## When the subroutine is executed
my %read_group = get_read_group(
    {
        infile_prefix    => $infile_prefix,
        platform         => $active_parameter{platform},
        sample_id        => q{ADM1059A1},
        sample_info_href => \%sample_info,
    }
);

## Then return expected read group hash
my %expected_read_group = (
    id   => q{ADM1059A1_161011_TestFilev2_GAGATTCC_lane1},
    lane => 1,
    lb   => q{ADM1059A1},
    pl   => q{ILLUMINA},
    pu   => q{TestFilev2.1.GAGATTC},
    sm   => q{ADM1059A1},
);
is_deeply( \%read_group, \%expected_read_group, q{Get read group} );

done_testing();
