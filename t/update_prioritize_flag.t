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
        q{MIP::Analysis}       => [qw{ update_prioritize_flag }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Analysis qw{ update_prioritize_flag };

diag(   q{Test update_prioritize_flag from Analysis.pm v}
      . $MIP::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my @recipes   = qw{ cnvnator_ar delly_call delly_reformat tiddit };
my %parameter = (
    cnvnator_ar    => { variant_caller => q{cnvnator}, },
    delly_call     => { variant_caller => q{delly}, },
    delly_reformat => { variant_caller => q{delly}, },
    manta          => { variant_caller => q{manta}, },
    tiddit         => { variant_caller => q{tiddit}, },
);
my $prioritize_key          = q{manta,delly,cnvnator,tiddit};
my $consensus_analysis_type = q{wes};

$prioritize_key = update_prioritize_flag(
    {
        consensus_analysis_type => $consensus_analysis_type,
        parameter_href          => \%parameter,
        prioritize_key          => $prioritize_key,
        recipes_ref             => \@recipes,
    }
);

is( $prioritize_key, q{manta}, q{Updated prioritize flag} );

done_testing();
