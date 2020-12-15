#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ set_gene_panel }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_gene_panel };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test set_gene_panel from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );
my $aggregate_gene_panel_file =
  catfile( $Bin, qw{ data 643594-miptest aggregated_gene_panel_test.txt } );
my $aggregate_gene_panels_key = q{select_file};
my $gene_panel                = q{TEST};
my $recipe_name_test          = q{vcfparser_ar};
my %sample_info;

my %header_info = (
    $gene_panel => {
        display_name => q{gene_panel_test},
        gene_panel   => $gene_panel,
        updated_at   => q{2016-12-08},
        version      => q{1.0},
    }
);

## Given no aggregate_gene_panel_file
set_gene_panel(
    {
        aggregate_gene_panels_key => $aggregate_gene_panels_key,
        recipe_name               => $recipe_name_test,
        sample_info_href          => \%sample_info,
    }
);

## Then don't set gene panel
is( $sample_info{$recipe_name_test}{$aggregate_gene_panels_key},
    undef, q{Leave gene panel unset if no gene panel file} );

## Given a bed like file with gene panel information
set_gene_panel(
    {
        aggregate_gene_panel_file => $aggregate_gene_panel_file,
        aggregate_gene_panels_key => $aggregate_gene_panels_key,
        recipe_name               => $recipe_name_test,
        sample_info_href          => \%sample_info,
    }
);

## Then set gene panel
is_deeply( $sample_info{$recipe_name_test}{$aggregate_gene_panels_key}{gene_panel},
    \%header_info, q{Set gene panel} );

## Given a not valid gene panel
trap {
    set_gene_panel(
        {
            aggregate_gene_panel_file =>
              catfile( $Bin, qw{ data test_data not_valid_gene_panel.bed } ),
            aggregate_gene_panels_key => $aggregate_gene_panels_key,
            recipe_name               => $recipe_name_test,
            sample_info_href          => \%sample_info,
        }
    )
};

## Then warn
like(
    $trap->stderr,
    qr/Unable \s+ to \s+ write \s+ select_file/xms,
    q{Throw warning log message}
);

## Given a file without gene_panel information
trap {
    set_gene_panel(
        {
            aggregate_gene_panel_file => catfile(
                $Bin,
                qw{ data references grch37_agilent_sureselect_targets_cre_-v1-.bed }
            ),
            aggregate_gene_panels_key => $aggregate_gene_panels_key,
            recipe_name               => $recipe_name_test,
            sample_info_href          => \%sample_info,
        }
    )
};

## Then throw fatal message and exit
like( $trap->stderr, qr/Unable \s to \s parse /xms, q{Throw fatal log message} );
like( $trap->die,    qr/Unable/xms,                 q{Die on failing regexp} );

done_testing();
