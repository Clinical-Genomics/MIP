#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
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
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Store}          => [qw{ parse_store_files }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Store qw{ parse_store_files };

diag(   q{Test parse_store_files from Store.pm v}
      . $MIP::Store::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given entries with duplicate paths
my @files = (
    {
        format     => q{bam},
        id         => q{id_1},
        path       => q{bamfile.bam},
        path_index => q{bamfile.bai},
        step       => q{gatk_baserecalibration},
        tag        => q{alignment_file},
    },
    {
        format     => q{bam},
        id         => q{id_2},
        path       => q{bamfile.bam},
        path_index => q{bamfile.bai},
        step       => q{gatk_baserecalibration},
        tag        => q{alignment_file},
    },
);
my $store_files_ref = parse_store_files(
    {
        store_files_ref => \@files,
    }
);

## Then remove the first one
my @expected_files = (
    {
        format     => q{bam},
        id         => q{id_2},
        path       => q{bamfile.bam},
        path_index => q{bamfile.bai},
        step       => q{gatk_baserecalibration},
        tag        => q{alignment_file},
    }
);
is_deeply( $store_files_ref, \@expected_files, q{Remove duplicate entries} );

done_testing();
