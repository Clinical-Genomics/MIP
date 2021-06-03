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
        q{MIP::Parameter}      => [qw{ get_order_of_parameters }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ get_order_of_parameters };

diag(   q{Test get_order_of_parameters from Parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a definition file
my @define_parameters_file_paths =
  ( catfile( $Bin, qw{ data test_data define_parameters.yaml } ) );

my @order_parameters = get_order_of_parameters(
    { define_parameters_files_ref => \@define_parameters_file_paths, } );

my @expected_order_parameters =
  qw{ mip cluster_constant_path case_id email_types bwa_mem bwa_mem_rapid_db load_env supported_capture_kit };

## Then the order of parameters as they appear in the definiton file should be added to the array
is_deeply( \@order_parameters, \@expected_order_parameters,
    q{Added order of parameters from definition file} );

done_testing();
