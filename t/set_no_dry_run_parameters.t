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
        q{MIP::Sample_info}    => [qw{ set_no_dry_run_parameters }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_no_dry_run_parameters };

diag(   q{Test set_no_dry_run_parameters from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );


## Given dry run
Readonly my $ANALYSIS_DATE => 1999;
Readonly my $MIP_VERSION => 1;
my %active_parameter = ( dry_run_all => 1, );

my %dry_run_info;

my %no_dry_run_info = (
    analysisrunstatus => q{not_finished},
    analysis_date     => $ANALYSIS_DATE,
    mip_version       => $MIP_VERSION,
);

my %sample_info;

set_no_dry_run_parameters(
    {
        analysis_date    => $ANALYSIS_DATE,
        is_dry_run_all   => $active_parameter{dry_run_all},
        mip_version      => $MIP_VERSION,
        sample_info_href => \%sample_info,
    }
);

## Then set no active parameters
is_deeply( \%sample_info, \%dry_run_info, q{Set no parameters} );

## Given no dry run
$active_parameter{dry_run_all} = 0;

set_no_dry_run_parameters(
    {
        analysis_date    => $ANALYSIS_DATE,
        is_dry_run_all   => $active_parameter{dry_run_all},
        mip_version      => $MIP_VERSION,
        sample_info_href => \%sample_info,
    }
);

## Then set true analysis parameters
is_deeply( \%sample_info, \%no_dry_run_info, q{Set active run parameters} );

done_testing();
