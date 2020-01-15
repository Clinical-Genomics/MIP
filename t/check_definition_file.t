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
        q{MIP::Definition}     => [qw{ check_definition_file }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Definition qw{ check_definition_file };

diag(   q{Test check_definition_file from Definition.pm v}
      . $MIP::Definition::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given definition files
my $level = q{mip};

my $definition_file_path =
  catfile( dirname($Bin), qw{ definitions mip_parameters.yaml } );

## Not required parameter definition keys to check
my $not_required_definition_file_path =
  catfile( dirname($Bin), qw{ definitions not_required_parameters.yaml } );

## required parameter definition keys to check
my $required_definition_file_path =
  catfile( dirname($Bin), qw{ definitions required_parameters.yaml } );

my %parameter = check_definition_file(
    {
        define_parameters_path            => $definition_file_path,
        not_required_definition_file_path => $not_required_definition_file_path,
        required_definition_file_path     => $required_definition_file_path,
    }
);

## Then a parameter hash should be returned
ok( exists $parameter{mip}, q{Checked parameter according to definition files} );

done_testing();
