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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log };

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Rankvariant} => [
            qw{ analysis_rankvariant analysis_rankvariant_unaffected analysis_rankvariant_sv analysis_rankvariant_sv_unaffected }
        ],
        q{MIP::Set::Analysis}  => [qw{ set_rankvariants_ar }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Analysis qw{ set_rankvariants_ar };
use MIP::Recipes::Analysis::Rankvariant
  qw{ analysis_rankvariant analysis_rankvariant_unaffected analysis_rankvariant_sv analysis_rankvariant_sv_unaffected };

diag(   q{Test set_rankvariants_ar from Analysis.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given no recipes and only unaffected
my $log                = test_log( {} );
my $sample_id          = q{sample_1};
my $sample_id_affected = q{sample_2};
my %active_parameter   = ( sample_ids => [$sample_id], );
my %parameter          = ( cache => { unaffected => [$sample_id], }, );
my %analysis_recipe;
my %expected_recipe = (
    rankvariant    => \&analysis_rankvariant_unaffected,
    sv_rankvariant => \&analysis_rankvariant_sv_unaffected,
);

## Set which recipe to use depending on sample affected status
set_rankvariants_ar(
    {
        analysis_recipe_href => \%analysis_recipe,
        log                  => $log,
        parameter_href       => \%parameter,
        sample_ids_ref       => $active_parameter{sample_ids},
    }
);

## Then set rankvariants unaffected recipes
is_deeply( \%analysis_recipe, \%expected_recipe,
    q{Set unaffected recipes when all samples are unaffected} );

## Given affected samples
push @{ $active_parameter{sample_ids} }, $sample_id_affected;

my %expected_affected_recipe = (
    rankvariant    => \&analysis_rankvariant,
    sv_rankvariant => \&analysis_rankvariant_sv,
);

## Set which recipe to use depending on sample affected status
set_rankvariants_ar(
    {
        analysis_recipe_href => \%analysis_recipe,
        log                  => $log,
        parameter_href       => \%parameter,
        sample_ids_ref       => $active_parameter{sample_ids},
    }
);

## Then set rankvaraints affected recipes
is_deeply( \%analysis_recipe, \%expected_affected_recipe,
    q{Set recipes when a sample has affected status} );

done_testing();
