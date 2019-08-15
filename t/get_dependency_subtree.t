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
        q{MIP::Get::Analysis}  => [qw{ get_dependency_subtree }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Analysis qw{ get_dependency_subtree };
use MIP::Test::Fixtures qw{ test_mip_hashes };

diag(   q{Test get_dependency_subtree from Analysis.pm v}
      . $MIP::Get::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %dependency_tree = test_mip_hashes(
    {
        mip_hash_name => q{dependency_tree_rna},
    }
);
my $chain_initiation_point  = q{CHAIN_QC};
my $dependency_subtree_href = {};

## Given request to get DELLY_CALL subtree
get_dependency_subtree(
    {
        chain_initiation_point  => $chain_initiation_point,
        dependency_tree_href    => \%dependency_tree,
        dependency_subtree_href => $dependency_subtree_href,
    },
);

## Then get it
my %expected_tree =
  ( CHAIN_QC => [ { PARALLEL => [qw{ preseq_ar rseqc genebody_coverage }] }, ] );
is_deeply( $dependency_subtree_href, \%expected_tree, q{Get subtree} );

done_testing();
