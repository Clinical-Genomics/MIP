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
        q{MIP::Sample_info}    => [qw{ get_sequence_run_type_is_interleaved }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ get_sequence_run_type_is_interleaved };

diag(   q{Test get_sequence_run_type_is_interleaved from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a not interleaved fastq file
my $infile_prefix = q{an_infile_prefix};
my $sample_id     = q{ADM1059A1};
my %sample_info   = (
    sample => {
        $sample_id => {
            file => {
                $infile_prefix => { sequence_run_type => { interleaved => undef, }, },
            },
        },
    },
);
my $is_interleaved = get_sequence_run_type_is_interleaved(
    {
        infile_lane_prefix => $infile_prefix,
        sample_id          => $sample_id,
        sample_info_href   => \%sample_info,
    }
);

## Then return false
is( $is_interleaved, undef, q{Got a not interleaved file} );

## Given an interleaved file
$sample_info{sample}{$sample_id}{file}{$infile_prefix}{sequence_run_type}{interleaved} =
  1;

$is_interleaved = get_sequence_run_type_is_interleaved(
    {
        infile_lane_prefix => $infile_prefix,
        sample_id          => $sample_id,
        sample_info_href   => \%sample_info,
    }
);

## Then return true
ok( $is_interleaved, q{Got an interleaved file} );

done_testing();
