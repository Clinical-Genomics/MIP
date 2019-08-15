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
our $VERSION = 1.01;

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
        q{MIP::Sample_info}    => [qw{ set_parameter_in_sample_info }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_parameter_in_sample_info };

diag(   q{Test set_parameter_in_sample_info from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
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
foreach my $key_to_add (@add_keys) {

    set_parameter_in_sample_info(
        {
            active_parameter_href => \%active_parameter,
            key_to_add            => $key_to_add,
            sample_info_href      => \%sample_info,
        }
    );
}

## Then add these to sample_info hash
is_deeply( \%active_parameter, \%sample_info, q{Added keys to sample info} );

done_testing();
