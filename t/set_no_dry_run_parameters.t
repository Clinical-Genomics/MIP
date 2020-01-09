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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = '1.0.0';

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
        q{MIP::Set::Parameter} => [qw{ set_no_dry_run_parameters }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Parameter qw{ set_no_dry_run_parameters };

diag(   q{Test set_no_dry_run_parameters from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given dry run
my $analysis_date = q{1999};

my %active_parameter = ( dry_run_all => 1, );

my %dry_run_info;

my %no_dry_run_info = (
    analysisrunstatus => q{not_finished},
    analysis_date     => $analysis_date,
    mip_version       => $VERSION,
);

my %sample_info;

set_no_dry_run_parameters(
    {
        is_dry_run_all   => $active_parameter{dry_run_all},
        analysis_date    => $analysis_date,
        mip_version      => $VERSION,
        sample_info_href => \%sample_info,
    }
);

## Then set no active parameters
is_deeply( \%sample_info, \%dry_run_info, q{Set no parameters} );

## Given no dry run
$active_parameter{dry_run_all} = 0;

set_no_dry_run_parameters(
    {
        is_dry_run_all   => $active_parameter{dry_run_all},
        analysis_date    => $analysis_date,
        mip_version      => $VERSION,
        sample_info_href => \%sample_info,
    }
);

## Then set true analysis parameters
is_deeply( \%sample_info, \%no_dry_run_info, q{Set active run parameters} );

done_testing();
