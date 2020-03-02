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
        q{MIP::Contigs}        => [qw{ set_contigs }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Contigs qw{ set_contigs };

diag(   q{Test set_contigs from List.pm v}
      . $MIP::Contigs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $GENOME_BUILD_38 => 38;
Readonly my $GENOME_BUILD_37 => 37;

my @grch38_contigs = qw{
  chr1 chr2 chr3 chr4 chr5 chr6
  chr7 chr8 chr9 chr10 chr11 chr12
  chr13 chr14 chr15 chr16 chr17 chr18
  chr19 chr20 chr21 chr22 chrX chrY
  chrM };
my @grch38_contigs_size_ordered = qw{
  chr1 chr2 chr3 chr4 chr5 chr6
  chr7 chrX chr8 chr9 chr10 chr11
  chr12 chr13 chr14 chr15 chr16 chr17
  chr18 chr19 chr20 chr21 chr22 chrY
  chrM };

my @grch37_contigs = qw{
  1 2 3 4 5 6 7 8 9 10
  11 12 13 14 15 16 17 18 19 20
  21 22 X Y MT };
my @grch37_contigs_size_ordered = qw{
  1 2 3 4 5 6 7 X 8 9
  10 11 12 13 14 15 16 17 18 19
  20 21 22 Y MT };

## Given a human genome reference source when 37
my %file_info = ( human_genome_reference_version => $GENOME_BUILD_37, );

set_contigs(
    {
        file_info_href => \%file_info,
        version        => $file_info{human_genome_reference_version},
    }
);

## Then contigs sets should be set for 37
is_deeply( \@{ $file_info{contigs} }, \@grch37_contigs, q{Set grch37 reference contigs} );
is_deeply(
    \@{ $file_info{contigs_size_ordered} },
    \@grch37_contigs_size_ordered,
    q{Set grch37 reference size ordered contigs}
);
is_deeply( \@{ $file_info{bam_contigs} },
    \@grch37_contigs, q{Set grch37 reference bam contigs} );
is_deeply(
    \@{ $file_info{bam_contigs_size_ordered} },
    \@grch37_contigs_size_ordered,
    q{Set grch37 reference size ordered bam contigs}
);

# Given a human genome reference source when 38 and alternative loci
$file_info{human_genome_reference_version} = $GENOME_BUILD_38;
my @grch38_alt_loci = qw{ chr1_KI270706v1_random
  chr1_KI270707v1_random };
$file_info{dict_contigs} = [ qw{ chr1 }, @grch38_alt_loci ];

set_contigs(
    {
        file_info_href => \%file_info,
        version        => $file_info{human_genome_reference_version},
    }
);

## Then contigs sets should be set for 38
is_deeply( \@{ $file_info{contigs} }, \@grch38_contigs, q{Set grch38 reference contigs} );
is_deeply(
    \@{ $file_info{contigs_size_ordered} },
    \@grch38_contigs_size_ordered,
    q{Set grch38 reference size ordered contigs}
);
is_deeply( \@{ $file_info{alt_loci} },
    \@grch38_alt_loci, q{Set grch38 reference alt loci} );

done_testing();
