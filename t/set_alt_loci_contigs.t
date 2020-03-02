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
        q{MIP::File_info}      => [qw{ set_alt_loci_contigs }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_alt_loci_contigs };

diag(   q{Test set_alt_loci_contigs from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given primary contigs and alternative loci
my $contig_set     = q{alt_loci};
my %primary_contig = (
    chr1 => undef,
    chr2 => undef,
);
my @grch38_alt_loci = qw{ chr1_KI270706v1_random
  chr1_KI270707v1_random };
my %file_info = ( dict_contigs => [ qw{ chr1 chr2 }, @grch38_alt_loci ], );

set_alt_loci_contigs(
    {
        alt_contig_set_name => $contig_set,
        file_info_href      => \%file_info,
        primary_contig_href => \%primary_contig,
    }
);

## Then
is_deeply( \@{ $file_info{$contig_set} },
    \@grch38_alt_loci, q{Set grch38 reference alt loci} );

done_testing();
