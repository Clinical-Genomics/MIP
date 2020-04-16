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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Parse::Gender}  => [qw{ parse_fastq_for_gender }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Gender qw{ parse_fastq_for_gender };

diag(   q{Test parse_fastq_for_gender from Gender.pm v}
      . $MIP::Parse::Gender::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1 } );

## Given no other gender
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => q{bwa_mem},
    }
);
my $consensus_analysis_type = q{wgs};
my $sample_id               = $active_parameter{sample_ids}[2];
$active_parameter{found_other} = 0;
my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
    }
);
push @{ $file_info{$sample_id}{mip_infiles} }, q{ADM1059A3.fastq};
my $infile_lane_prefix = $sample_id;
my %infile_lane_prefix = ( $sample_id => [ $infile_lane_prefix, $infile_lane_prefix ], );

my %sample_info = (
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

my $is_gender_other = parse_fastq_for_gender(
    {
        active_parameter_href   => \%active_parameter,
        consensus_analysis_type => $consensus_analysis_type,
        file_info_href          => \%file_info,
        infile_lane_prefix_href => \%infile_lane_prefix,
        sample_info_href        => \%sample_info,
    }
);

## Then skip estimation using reads
is( $is_gender_other, undef, q{No unknown gender} );

## Given a sample when gender unknown
$active_parameter{found_other}    = 1;
$active_parameter{found_male}     = 1;
$active_parameter{gender}{others} = [$sample_id];

$is_gender_other = parse_fastq_for_gender(
    {
        active_parameter_href   => \%active_parameter,
        consensus_analysis_type => $consensus_analysis_type,
        file_info_href          => \%file_info,
        infile_lane_prefix_href => \%infile_lane_prefix,
        sample_info_href        => \%sample_info,
    }
);

## Then skip estimation using reads
is( $is_gender_other, undef, q{Unknown gender} );

done_testing();
