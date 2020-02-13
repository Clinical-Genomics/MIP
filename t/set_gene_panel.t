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
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ set_gene_panel }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_gene_panel };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test set_gene_panel from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
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
my $case_id_test              = q{case_id};
my $recipe_name_test          = q{vcfparser_ar};
my %sample_info;

my %header_info = (
    display_name => q{gene_panel_test},
    gene_panel   => $gene_panel,
    updated_at   => q{2016-12-08},
    version      => q{1.0},
);

set_gene_panel(
    {
        aggregate_gene_panel_file => $aggregate_gene_panel_file,
        aggregate_gene_panels_key => $aggregate_gene_panels_key,
        case_id                   => $case_id_test,
        log                       => $log,
        recipe_name               => $recipe_name_test,
        sample_info_href          => \%sample_info,
    }
);

is(
    exists $sample_info{$recipe_name_test}{$aggregate_gene_panels_key}
      {gene_panel}{$gene_panel},
    1,
    q{Gene panel key added to $sample_info}
);

while ( my ( $key, $value ) = each %header_info ) {

## Test gene panel info
    my $set_header_value =
      $sample_info{$recipe_name_test}{$aggregate_gene_panels_key}{gene_panel}
      {$gene_panel}{$key};

    is( $set_header_value, $value,
        q{Gene panel header info value for key: } . $key . q{ added to $sample_info} );
}

## Given a not valid gene panel
trap {
    set_gene_panel(
        {
            aggregate_gene_panel_file =>
              catfile( $Bin, qw{ data test_data not_valid_gene_panel.bed } ),
            aggregate_gene_panels_key => $aggregate_gene_panels_key,
            case_id                   => $case_id_test,
            log                       => $log,
            recipe_name               => $recipe_name_test,
            sample_info_href          => \%sample_info,
        }
    )
};

## Then warn
like(
    $trap->stderr,
    qr/Unable \s+ to \s+ write \s+ select_file/xms,
    q{Throw fatal log message}
);

done_testing();
