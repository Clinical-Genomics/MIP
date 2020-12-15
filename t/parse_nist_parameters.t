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
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Reference}      => [qw{ parse_nist_parameters }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ parse_nist_parameters };

diag(   q{Test parse_nist_parameters from Reference.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given nist info
my %active_parameter = (
    nist_call_set_vcf =>
      { q{3.3.2} => { NA12878 => q{grch37_nist_hg001_-na12878_v3.3.2-.vcf}, }, },
    nist_call_set_bed =>
      { q{3.3.2} => { NA12878 => q{grch37_nist_hg001_-na12878_v3.3.2-.bed}, }, },
    nist_id       => { sample_1 => q{NA12878}, },
    nist_versions => [qw{ 3.3.2 }],
    reference_dir => catdir( $Bin, qw{ data references } ),
    rtg_vcfeval   => 1,
    sample_ids    => [qw{ sample_1 }],
);

my $is_ok = parse_nist_parameters(
    {
        active_parameter_href => \%active_parameter,
    }
);

## Then return true
ok( $is_ok, q{Parsed nist parameters} );

done_testing();
