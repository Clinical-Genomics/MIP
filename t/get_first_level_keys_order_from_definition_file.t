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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Definition}     => [qw{ get_first_level_keys_order_from_definition_file }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Definition qw{ get_first_level_keys_order_from_definition_file };

diag(   q{Test get_first_level_keys_order_from_definition_file from Definition.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $yaml_file = catfile( $Bin, qw{ data test_data define_parameters.yaml } );

my @order_parameters = get_first_level_keys_order_from_definition_file(
    {
        file_path => $yaml_file
    }
);

my $first_key             = $order_parameters[0];
my $cluster_constant_path = $order_parameters[1];

isnt( $first_key, q{---}, q{Skipped header and comments} );

is( $first_key, q{mip}, q{Found first key} );

is( $cluster_constant_path, q{cluster_constant_path}, q{Found second key} );

done_testing();
