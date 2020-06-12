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
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Download}       => [qw{ get_download_reference_attributes }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Download qw{ get_download_reference_attributes };

diag(   q{Test get_download_reference_attributes from Download.pm v}
      . $MIP::Download::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given some references
my %active_parameter =
  test_mip_hashes( { mip_hash_name => q{download_active_parameter}, } );
my $reference_id   = q{human_reference};
my $genome_version = q{grch38};
my $version        = q{assembly38};

## When getting a reference version
my %reference = get_download_reference_attributes(
    {
        id             => $reference_id,
        genome_version => $genome_version,
        reference_href => $active_parameter{reference_feature},
        version        => $version,
    }
);

## Then reference file hash should be returned
is_deeply(
    \%reference,
    $active_parameter{reference_feature}{$reference_id}{$genome_version}{$version},
    q{Got reference version attributes}
);

## When not supplying a version
my $return = get_download_reference_attributes(
    {
        id             => $reference_id,
        genome_version => $genome_version,
        reference_href => $active_parameter{reference_feature},
        version        => q{does not exists},
    }
);
is( $return, undef, q{Did not get reference version attributes} );

done_testing();
