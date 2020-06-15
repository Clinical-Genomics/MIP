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
        q{MIP::Vcfanno}        => [qw{ get_toml_annotation_preops }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfanno qw{ get_toml_annotation_preops };

diag(   q{Test get_toml_annotation_preops from Vcfanno.pm v}
      . $MIP::Vcfanno::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given no preops
my %annotation = (
    file   => q{reference.vcf},
    fields => [q{vcf_tag}],
    ops    => [q{self}],
    names  => [q{vcf_tag}],
);

my %preops = get_toml_annotation_preops(
    {
        annotation_href => \%annotation,
    }
);

## Then return undef
is( %preops, 0, q{Return empty hash on missing preops} );

## Given preops
$annotation{preops} =
  { fields => [q{vcf_tag}], ops => [q{some_code.lua}], names => [q{mod_tag}], };

%preops = get_toml_annotation_preops(
    {
        annotation_href => \%annotation,
    }
);

## Then return preops annotation
my %expected = (
    file   => q{reference.vcf},
    fields => [q{vcf_tag}],
    ops    => [q{some_code.lua}],
    names  => [q{mod_tag}],
);
is_deeply( \%preops, \%expected, q{Return preops} );

done_testing();
