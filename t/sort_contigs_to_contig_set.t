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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Contigs} => [qw{ sort_contigs_to_contig_set }],
        q{MIP::Test::Fixtures}  => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Contigs qw{ sort_contigs_to_contig_set };

diag(   q{Test sort_contigs_to_contig_set from Contigs.pm v}
      . $MIP::Contigs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

my $consensus_analysis_type = q{wes};

## Given a reference hash of array, when no hash of array to sort exists
my %file_info = ( contigs_size_ordered => [qw{ chr1 chr2 chr3 chrM}],
select_file_contigs => [],
);

trap {
    sort_contigs_to_contig_set(
        {
            consensus_analysis_type => $consensus_analysis_type,
            sort_contigs_ref        => $file_info{select_file_contigs},
            sort_reference_contigs_ref => $file_info{contigs_size_ordered},
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if hash key to sort does not exist in supplied hash} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if hash key to sort does not exist in supplied hash} );

## Given a hash of array to sort, when no hash of array as reference exists
%file_info = ( contigs_size_ordered => [],
select_file_contigs => [qw{chr2 chr1 chrM chr3}], );

trap {
    sort_contigs_to_contig_set(
        {
            consensus_analysis_type => $consensus_analysis_type,
            sort_contigs_ref        => $file_info{select_file_contigs},
            sort_reference_contigs_ref => $file_info{contigs_size_ordered},
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if hash key for reference does not exist in supplied hash} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if hash key for reference does not exist in supplied hash}
);

## Given wes, when lacking mitochondrial contig in array to sort
%file_info = (
    contigs_size_ordered => [qw{ chr1 chr2 chr3 chrM}],
    select_file_contigs  => [qw{chr2 chr1 chr3}],
);

@{ $file_info{sorted_select_file_contigs} } = sort_contigs_to_contig_set(
    {
        consensus_analysis_type => $consensus_analysis_type,
        sort_contigs_ref        => $file_info{select_file_contigs},
        sort_reference_contigs_ref => $file_info{contigs_size_ordered},
    }
);

## Remove chrM for comparison
pop @{ $file_info{contigs_size_ordered} };

## Then return sorted array according to reference (minus chrM in ref)
is_deeply(
    \@{ $file_info{contigs_size_ordered} },
    \@{ $file_info{sorted_select_file_contigs} },
    q{Sorted contigs according to reference for wes}
);

## Given wes, when lacking a contig in reference array
%file_info = (
    contigs_size_ordered => [qw{ chr2 chr3 chrM}],
    select_file_contigs  => [qw{chr2 chr1 chrM chr3}],
);

trap {
    sort_contigs_to_contig_set(
        {
            consensus_analysis_type => $consensus_analysis_type,
            sort_contigs_ref        => $file_info{select_file_contigs},
            sort_reference_contigs_ref => $file_info{contigs_size_ordered},
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if contig is lacking from sorted array due to reference} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if contig is lacking from sorted array due to reference} );

## Given all ok parameters
%file_info = (
    contigs_size_ordered => [qw{ chr1 chr2 chr3 chrM}],
    select_file_contigs  => [qw{chr2 chr1 chrM chr3}],
);

@{ $file_info{sorted_select_file_contigs} } = sort_contigs_to_contig_set(
    {
        consensus_analysis_type => $consensus_analysis_type,
        sort_contigs_ref        => $file_info{select_file_contigs},
        sort_reference_contigs_ref => $file_info{contigs_size_ordered},
    }
);

## Then return a complete sorted array
is_deeply(
    \@{ $file_info{contigs_size_ordered} },
    \@{ $file_info{sorted_select_file_contigs} },
    q{Sorted contigs according to reference}
);

done_testing();
