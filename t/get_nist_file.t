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
        q{MIP::Reference}      => [qw{ get_nist_file }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ get_nist_file };

diag(   q{Test get_nist_file from Reference.pm v}
      . $MIP::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given nist info
my %nist = ( nist_call_set_vcf =>
      { q{3.3.2} => { NA12878 => q{grch37_nist_hg001_-na12878_v3.3.2-.vcf}, }, }, );
my $nist_id        = q{NA12878};
my $nist_parameter = q{nist_call_set_vcf};
my $nist_version   = q{3.3.2};

my $nist_file = get_nist_file(
    {
        nist_href    => $nist{$nist_parameter},
        nist_id      => $nist_id,
        nist_version => $nist_version,
    }
);

my $expected_file = q{grch37_nist_hg001_-na12878_v3.3.2-.vcf};

## Then nist file should be returned
is( $nist_file, $expected_file, q{Got nist file} );

done_testing();
