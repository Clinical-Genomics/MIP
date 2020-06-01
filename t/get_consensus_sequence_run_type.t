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
        q{MIP::File_info}      => [qw{ get_consensus_sequence_run_type }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ get_consensus_sequence_run_type };

diag(   q{Test get_consensus_sequence_run_type from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given file info
my $mip_file_format = q{ADM1059A1_161011_HHJJCCCXY_NAATGCGC_lane1};
my $sample_id       = q{ADM1059A1};
my @sample_ids      = ( $sample_id, );
my %file_info       = test_mip_hashes( { mip_hash_name => q{file_info}, } );

## When all files have the same sequence_run_type
my $consensus_type = get_consensus_sequence_run_type(
    {
        file_info_href => \%file_info,
        sample_ids_ref => \@sample_ids,
    }
);

## Then return consensus type
is( $consensus_type, q{single-end}, q{Returned consensus sequence run type} );

## When multiple sequence run types
push @{ $file_info{$sample_id}{no_direction_infile_prefixes} }, $mip_file_format;
$file_info{$sample_id}{$mip_file_format}{sequence_run_type} = q{paired-end};

my $has_consensus = get_consensus_sequence_run_type(
    {
        file_info_href => \%file_info,
        sample_ids_ref => \@sample_ids,
    }
);

## Then return false
is( $has_consensus, 0, q{Found different sequence run type} );

done_testing();
