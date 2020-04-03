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
        q{MIP::Active_parameter} => [qw{ set_nist_file_name_path }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_nist_file_name_path };

diag(   q{Test set_nist_file_name_path from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given nist info
my %active_parameter = (
    nist_call_set_vcf =>
      { q{3.3.2} => { NA12878 => q{grch37_nist_hg001_-na12878_v3.3.2-.vcf}, }, },
    nist_call_set_bed =>
      { q{3.3.2} => { NA12878 => q{grch37_nist_hg001_-na12878_v3.3.2-.bed}, }, },
    nist_id       => { sample_1 => q{NA12878}, },
    nist_versions => [qw{ 3.3.2 }],
    reference_dir => catdir( $Bin, qw{ data references } ),
    sample_ids    => [qw{ sample_1 }],
);
my @nist_parameters = (qw{ nist_call_set_vcf nist_call_set_bed });

my $is_ok = set_nist_file_name_path(
    {
        active_parameter_href => \%active_parameter,
        nist_parameters_ref   => \@nist_parameters,
    }
);

## Then return true
ok( $is_ok, q{Set nist file name path in nist hash parameters} );

## Then reference dir should have been prepended to file name
my $expected_path =
  catdir( $Bin, qw{ data references grch37_nist_hg001_-na12878_v3.3.2-.vcf} );

is( $active_parameter{nist_call_set_vcf}{q{3.3.2}}{NA12878},
    $expected_path, q{Set path using reference dir} );

done_testing();
