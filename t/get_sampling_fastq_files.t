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
        q{MIP::Parse::Gender}  => [qw{ get_sampling_fastq_files }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Gender qw{ get_sampling_fastq_files };

diag(   q{Test get_sampling_fastq_files from Gender.pm v}
      . $MIP::Parse::Gender::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => q{bwa_mem},
    }
);
my $sample_id = $active_parameter{sample_ids}[2];

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

my ( $is_interleaved_fastq, @fastq_files ) = get_sampling_fastq_files(
    {
        infile_lane_prefix_href => \%infile_lane_prefix,
        infile_paths_ref        => $file_info{$sample_id}{mip_infiles},
        sample_id               => $sample_id,
        sample_info_href        => \%sample_info,
    }
);
my @expected_fastq_files = qw{ ADM1059A3.fastq ADM1059A3.fastq };
## Then return undef for interleaved
is( $is_interleaved_fastq, undef, q{No interleaved files} );

## Then return fastq files
is_deeply( \@fastq_files, \@expected_fastq_files, q{Got fastq files} );

done_testing();
