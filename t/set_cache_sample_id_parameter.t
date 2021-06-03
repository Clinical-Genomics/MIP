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
        q{MIP::Parameter}      => [qw{ set_cache_sample_id_parameter }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ set_cache_sample_id_parameter };

diag(   q{Test set_cache_sample_id_parameter from Parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a parameter and sample_id
my %parameter;
my $parameter_name  = q{plink_phenotype};
my $parameter_value = q{other};
my $sample_id       = q{sample_1};

set_cache_sample_id_parameter(
    {
        parameter_href  => \%parameter,
        parameter_name  => $parameter_name,
        parameter_value => $parameter_value,
        sample_id       => $sample_id,
    }
);

## Then set cache at sample level for parameter
is( $parameter{cache}{$sample_id}{$parameter_name},
    $parameter_value, q{Set parameter at sample level in cache} );

done_testing();
