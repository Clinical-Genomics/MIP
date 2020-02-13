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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ get_sequence_run_type }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ get_sequence_run_type };

diag(   q{Test get_sequence_run_type from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Define test hashes
my $sample_id               = q{test};
my $infile_lane_prefix      = q{test_lane1};
my %infile_lane_prefix_hash = ( $sample_id => [$infile_lane_prefix], );
my %sample_info             = (
    sample => {
        $sample_id => {
            file => {
                $infile_lane_prefix => {
                    sequence_run_type => q{paired-end},
                },
            },
        },
    },
);

## Given the test hashes and a sample id
my %sequence_run_type = get_sequence_run_type(
    {
        infile_lane_prefix_href => \%infile_lane_prefix_hash,
        sample_id               => $sample_id,
        sample_info_href        => \%sample_info,
    }
);

## Then return a hash with sequence run type per lane
my %expected = ( $infile_lane_prefix => q{paired-end} );
is_deeply( \%sequence_run_type, \%expected, q{Return sequence run type hash} );

## Given sample_info test hash, infile_prefix and sample_id
my $sequence_run_type_value = get_sequence_run_type(
    {
        infile_lane_prefix => $infile_lane_prefix,
        sample_id          => $sample_id,
        sample_info_href   => \%sample_info,
    }
);

## Then return the sequence run type for that lane
is( $sequence_run_type_value, q{paired-end}, q{Return sequence run type} );

## Given missing infile_lane_prefix information
trap {
    $sequence_run_type_value = get_sequence_run_type(
        {
            sample_id        => $sample_id,
            sample_info_href => \%sample_info,
        }
    )
};
## Then croak
is( $trap->leaveby, q{die}, q{Exit if lane info is missing} );

done_testing();
