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
        q{MIP::Vcfanno}        => [qw{ check_toml_annotation_for_tags }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfanno qw{ check_toml_annotation_for_tags };

diag(   q{Test check_toml_annotation_for_tags from Vcfanno.pm v}
      . $MIP::Vcfanno::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $bcftools_binary_path = q{bcftools};
my %missing_tag;
my %annotation = (
    file => catfile(
        $Bin, qw{ data references grch37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz }
    ),
    fields => [q{EUR_AF}],
    ops    => [q{self}],
    names  => [q{EUR_AF}],
);

my $is_ok = check_toml_annotation_for_tags(
    {
        annotation_href      => \%annotation,
        bcftools_binary_path => $bcftools_binary_path,
        missing_tag_href     => \%missing_tag,
    }
);

## Then return ok;
ok( $is_ok, q{Names and fields match} );

## Given missing tag
$annotation{fields} = [q{MISSING_TAG}];
$annotation{names}  = [q{MISSING_TAG}];
%missing_tag        = ();

check_toml_annotation_for_tags(
    {
        annotation_href      => \%annotation,
        bcftools_binary_path => $bcftools_binary_path,
        missing_tag_href     => \%missing_tag,
    }
);

## Then populate missing tag hash
my %expected = (
    catfile( $Bin,
        qw{ data references grch37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz } ) =>
      [qw{ MISSING_TAG}], );

is_deeply( \%missing_tag, \%expected, q{Populate missing tag hash} );

done_testing();
