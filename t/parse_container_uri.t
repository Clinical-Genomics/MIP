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
    my %perl_module = ( q{MIP::Environment::Container} => [qw{ parse_container_uri }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Container qw{ parse_container_uri };

diag(   q{Test parse_container_uri from Container.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a dockerhub uri
my $uri = q{docker.io/clinicalgenomics/chanjo:4.2.0};
## When container manager is singularity
parse_container_uri(
    {
        container_manager => q{singularity},
        uri_ref           => \$uri,
    }
);

## Then prepend docker://
my $expected_uri = q{docker://docker.io/clinicalgenomics/chanjo:4.2.0};
is( $uri, $expected_uri, q{Parse uri for singularity} );

## Given a dockerhub uri
$uri = q{docker.io/clinicalgenomics/chanjo:4.2.0};
## When container manager is docker
parse_container_uri(
    {
        container_manager => q{docker},
        uri_ref           => \$uri,
    }
);

## Then leave uri unchanged
$expected_uri = q{docker.io/clinicalgenomics/chanjo:4.2.0};
is( $uri, $expected_uri, q{Parse uri for docker} );

## Given a quay uri
$uri = q{quay.io/clinicalgenomics/chanjo:4.2.0};

## When container manager is singularity
parse_container_uri(
    {
        container_manager => q{singularity},
        uri_ref           => \$uri,
    }
);

## Then leave uri unchanged
$expected_uri = q{docker://quay.io/clinicalgenomics/chanjo:4.2.0};
is( $uri, $expected_uri, q{Parse quay uri for singularity} );

## Given a ghcr uri
my $uri = q{ghcr.io/dnanexus-rnd/glnexus:v1.4.1};
## When container manager is singularity
parse_container_uri(
    {
        container_manager => q{singularity},
        uri_ref           => \$uri,
    }
);

## Then prepend docker://
my $expected_uri = q{docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1};
is( $uri, $expected_uri, q{Parse uri for singularity} );

done_testing();
