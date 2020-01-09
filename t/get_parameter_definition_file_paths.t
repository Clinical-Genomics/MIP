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
use MIP::Constants qw{ $COMMA $SPACE $UNDERSCORE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Definition}     => [qw{ get_parameter_definition_file_paths }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Definition qw{ get_parameter_definition_file_paths };

diag(   q{Test get_parameter_definition_file_paths from Definition.pm v}
      . $MIP::Definition::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a first level definition
my $level                          = q{mandatory};
my $mandatory_definition_file_path = get_parameter_definition_file_paths( { level => $level, } );

my $expected_mandatory_file_path =
  catfile( dirname($Bin), qw{ t definitions mandatory_parameters.yaml } );

## Then get the definition file paths
is( $mandatory_definition_file_path,
    $expected_mandatory_file_path, q{Got mandatory definition file} );

## Given a third level definition
$level = q{rd_dna};
my @rd_dna_definition_file_paths = get_parameter_definition_file_paths( { level => $level, } );

my @expected_rd_dna_file_paths = (
    catfile( dirname($Bin), qw{ t definitions mip_parameters.yaml } ),
    catfile( dirname($Bin), qw{ t definitions analyse_parameters.yaml } ),
    catfile( dirname($Bin), qw{ t definitions rd_dna_parameters.yaml } ),
);

## Then get the definition file paths
is_deeply(
    \@rd_dna_definition_file_paths,
    \@expected_rd_dna_file_paths,
    q{Got mip analyse rd_dna definition files}
);

done_testing();
