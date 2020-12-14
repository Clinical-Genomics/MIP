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
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Qccollect}      => [qw{ plink_relation_check }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qccollect qw{ plink_relation_check };

diag(   q{Test plink_relation_check from Qccollect.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample when not related
my $recipe_name = q{plink_relation_check_infile};
my $sample_id   = q{sample_1};
my %qc_data     = ( sample => { $sample_id => { $recipe_name => q{FAIL}, }, } );

plink_relation_check(
    {
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        sample_id    => $sample_id,
    }
);

my $status =
  q{Status:} . $recipe_name . q{:} . $qc_data{sample}{$sample_id}{$recipe_name};

## Then
is_deeply( \@{ $qc_data{evaluation}{$recipe_name} },
    [$status], q{Added FAILED evaluation status} );

done_testing();
