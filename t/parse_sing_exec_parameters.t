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
use Modern::Perl qw{ 2017 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Parse::Parameter} => [qw{ parse_sing_exec_parameters }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Parameter qw{ parse_sing_exec_parameters };

diag(   q{Test parse_sing_exec_parameter from Parse::Parameter.pm v}
      . $MIP::Parse::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given some input parameters were some are existing paths.
my @test_parameters = (
    catdir( dirname($Bin), qw{ t data modules miniconda } ),
    catfile( dirname($Bin), qw{ templates mip_rd_dna_config.yaml } ),
    catdir( dirname($Bin), qw{ t data modules mumbo jumbo } ),
    catdir(qw{ dir_a dir_c }),
    qw{ some_option an_executable . },
);

my @return_paths =
  parse_sing_exec_parameters( { sing_exec_parameters_ref => \@test_parameters, } );

## Then return existing paths and remove parts of paths that doesn't exist
my @expected = (
    catdir( dirname($Bin), qw{ t data modules miniconda } ),
    catfile( dirname($Bin), qw{ templates mip_rd_dna_config.yaml } ),
    catdir( dirname($Bin), qw{ t data modules } ),
);

is_deeply( \@return_paths, \@expected, q{Return paths} );

done_testing();
