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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ set_in_sample_info }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_in_sample_info };

diag(   q{Test set_in_sample_info from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my @add_keys = qw{ analysis_type expected_coverage };

my %active_parameter = (
    analysis_type => {
        sample_1 => q{wgs},
        sample_2 => q{wes},
    },
    expected_coverage => {
        sample_1 => 1,
        sample_2 => 1,
    },
);
my %sample_info;

## Given parameters in active_parameters
PARAMETER:
foreach my $key (@add_keys) {

    set_in_sample_info(
        {
            key              => $key,
            sample_info_href => \%sample_info,
            value            => \%{ $active_parameter{$key} },
        }
    );
}

## Then add these to sample_info hash
is_deeply( \%active_parameter, \%sample_info, q{Added keys to sample info} );

done_testing();
