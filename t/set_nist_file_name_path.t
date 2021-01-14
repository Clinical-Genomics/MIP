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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Reference}      => [qw{ set_nist_file_name_path }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ set_nist_file_name_path };

diag(   q{Test set_nist_file_name_path from Reference.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given nist info
my %nist = ( nist_call_set_vcf =>
      { q{3.3.2} => { NA12878 => q{grch37_nist_hg001_-na12878_v3.3.2-.vcf}, }, }, );
my $file_name      = q{grch37_nist_hg001_-na12878_v3.3.2-.vcf};
my $nist_id        = q{NA12878};
my $nist_parameter = q{nist_call_set_vcf};
my $nist_version   = q{3.3.2};
my $reference_dir  = catdir( $Bin, qw{ data references } );

my $is_ok = set_nist_file_name_path(
    {
        file_name     => $file_name,
        nist_href     => $nist{nist_call_set_vcf},
        nist_id       => $nist_id,
        nist_version  => $nist_version,
        reference_dir => $reference_dir,
    }
);

## Then return true
ok( $is_ok, q{Set nist file name path in nist hash parameters} );

## Then reference dir should have been prepended to file name
my $expected_path = catdir( $reference_dir, qw{ grch37_nist_hg001_-na12878_v3.3.2-.vcf} );

is( $nist{nist_call_set_vcf}{q{3.3.2}}{NA12878},
    $expected_path, q{Set path using reference dir} );

done_testing();
