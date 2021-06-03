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
        q{MIP::Active_parameter} => [qw{ set_pedigree_sample_id_parameter }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_pedigree_sample_id_parameter };

diag(   q{Test set_pedigree_sample_id_parameter from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a pedigree parameter
my %active_parameter;
my $pedigree_key   = q{phenotype};
my $pedigree_value = q{affected};
my $sample_id      = q{sample_1};

set_pedigree_sample_id_parameter(
    {
        active_parameter_href => \%active_parameter,
        pedigree_key          => $pedigree_key,
        pedigree_value        => $pedigree_value,
        sample_id             => $sample_id,
    }
);

## Then pedigree parameter should be set in active_parameter
is( $active_parameter{phenotype}{$sample_id}, q{affected}, q{Set pedigree parameter} );

done_testing();
