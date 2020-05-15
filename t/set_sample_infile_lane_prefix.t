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
        q{MIP::File_info}      => [qw{ set_sample_infile_lane_prefix }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_sample_infile_lane_prefix };

diag(   q{Test set_sample_infile_lane_prefix from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a lane_tracker, sample_id and an infile_lane_prefix
my %file_info;
my $lane_tracker       = 0;
my $infile_lane_prefix = q{infile_lane_prefix};
my $sample_id          = q{sample_id};

## When direction is not one
my $direction = 2;

## Then do not set infile_lane_prefix in file_info hash
set_sample_infile_lane_prefix(
    {
        direction       => $direction,
        file_info_href  => \%file_info,
        lane_tracker    => $lane_tracker,
        mip_file_format => $infile_lane_prefix,
        sample_id       => $sample_id,
    }
);

is( @{ $file_info{$sample_id}{infile_lane_prefix} },
    0, q{Did not set infile_prefix for read direction two} );

## When direction is one
$direction = 1;

## Then set infile_lane_prefix in file_info hash
set_sample_infile_lane_prefix(
    {
        direction       => $direction,
        file_info_href  => \%file_info,
        lane_tracker    => $lane_tracker,
        mip_file_format => $infile_lane_prefix,
        sample_id       => $sample_id,
    }
);

is( @{ $file_info{$sample_id}{infile_lane_prefix} },
    1, q{Set infile_prefix for read direction one} );

done_testing();
