#! /usr/bin/env perl

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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Contigs}        => [qw{ get_contig_set }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Contigs qw{ get_contig_set };
use MIP::Test::Fixtures qw { test_mip_hashes };

diag(   q{Test get_contig_set from Contigs.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $GENOME_VERSION => 38;

## Given a version that does not exists
my $faulty_version = q{this_is_not_a_genome_version};

## When getting contigs using the faulty version
my $result = get_contig_set( { version => $faulty_version, } );

## Then return undef
is_deeply( $result, undef, q{Version does not exists in contig hash} );

## Given a version

## When getting contigs using genome version
my %contigs = get_contig_set( { version => $GENOME_VERSION, } );

my %expected_contigs = test_mip_hashes( { mip_hash_name => q{primary_contig}, } );

## Then return all contig sets for version
is_deeply(
    \%contigs,
    $expected_contigs{$GENOME_VERSION},
    q{Got primary contig sets for genome version}
);

## Given a version and a set with a array ref
my $array_set = q{contigs};

## When getting contigs using genome version and a set
my @contigs = get_contig_set(
    {
        contig_set => $array_set,
        version    => $GENOME_VERSION,
    }
);

## Then return contig set
is_deeply(
    \@contigs,
    $expected_contigs{$GENOME_VERSION}{$array_set},
    q{Got primary contig array set for genome version and set}
);

## Given a version and a set with a hash ref
my $hash_set = q{synonyms_map};

## When getting contigs using genome version and a set
my %contigs_synonyms_map = get_contig_set(
    {
        contig_set => $hash_set,
        version    => $GENOME_VERSION,
    }
);

## Then return contig set
is_deeply(
    \%contigs_synonyms_map,
    $expected_contigs{$GENOME_VERSION}{$hash_set},
    q{Got primary contig hash set for genome version and set}
);

done_testing();
