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
use MIP::Constants
  qw{ $BACKWARD_SLASH $COLON $COMMA $DOUBLE_QUOTE $FORWARD_SLASH $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.04;

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
        q{MIP::Recipes::Install::Cadd} => [qw{ install_cadd }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Install::Cadd qw{ install_cadd };

diag(   q{Test install_cadd from Cadd.pm v}
      . $MIP::Recipes::Install::Cadd::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given install parameters
my %active_parameter =
  test_mip_hashes( { mip_hash_name => q{install_rd_dna_active_parameter}, } );
$active_parameter{reference_dir} = catdir(qw{ reference dir });
my $is_ok = install_cadd(
    {
        active_parameter_href => \%active_parameter,
        container_href        => $active_parameter{singularity}{cadd},
        container_path        => catfile(q{cadd.sif}),
        filehandle            => $filehandle,
    }
);

## Then return TRUE
ok( $is_ok, q{Executed install cadd recipe } );

## Given no reference directory
$active_parameter{reference_dir} = undef;
$active_parameter{singularity}{cadd}{program_bind_paths} = undef;
my $expected = [
    catdir(
        $BACKWARD_SLASH
          . $DOUBLE_QUOTE
          . $BACKWARD_SLASH
          . q{$MIP_BIND}
          . $BACKWARD_SLASH
          . $DOUBLE_QUOTE,
        qw{ CADD-scripts data annotations }
      )
      . $COLON
      . catdir( $FORWARD_SLASH, qw{ opt CADD-scripts data annotations } )
];
install_cadd(
    {
        active_parameter_href => \%active_parameter,
        container_href        => $active_parameter{singularity}{cadd},
        container_path        => catfile(q{cadd.sif}),
        filehandle            => $filehandle,
    }
);

## Then append relative bind baths
is_deeply( $active_parameter{singularity}{cadd}{program_bind_paths},
    $expected, q{Add MIP_BIND relative bind path} );
close $filehandle;

done_testing();
