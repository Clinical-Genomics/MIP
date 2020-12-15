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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE $UNDERSCORE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Definition}     => [qw{ get_dependency_tree_from_definition_file }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Definition qw{ get_dependency_tree_from_definition_file };

diag(   q{Test get_dependency_tree_from_definition_file from Definition.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a level definition
my $level = q{rd_dna};
my %dependency_tree = get_dependency_tree_from_definition_file( { level => $level, } );

my %expected_dependency_tree = (CHAIN_ALL => [ {CHAIN_FASTQ => [qw{ fastqc_ar }]},
],);

## Then get the dependency tree initiation map
is_deeply(
    \%dependency_tree,
    \%expected_dependency_tree,
    q{Got mip analyse rd_dna initiation map hash}
);


done_testing();
