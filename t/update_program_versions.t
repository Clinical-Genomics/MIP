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
use autodie qw{ :all };
use Clone qw{ clone };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Script::Utils}  => [qw{ update_program_versions }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Utils qw{ update_program_versions };

diag(   q{Test update_program_versions from Utils.pm v}
      . $MIP::Script::Utils::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given request for version updates
my $parameter_href = {
    installations    => [q{emip}],
    program_versions => {
        samtools => q{1.9},
        bedtools => q{2.30},
    },
    emip => {
        conda => {
            samtools => q{1.8},
        },
        shell => {
            bedtools => {
                version => q{2.27},
            },
        },
    },
};

update_program_versions(
    {
        parameter_href => $parameter_href,
    }
);

## Then update and remove program_versions key
my $expected_href = {
    installations => [q{emip}],
    emip          => {
        conda => {
            samtools => q{1.9},
        },
        shell => {
            bedtools => {
                version => q{2.30},
            },
        },
    },
};
is_deeply( $parameter_href, $expected_href, q{Update program versions} );

## Given no program versions to update
my $test_href = clone($expected_href);
update_program_versions(
    {
        parameter_href => $test_href,
    }
);

## Then return unchanged hash
is_deeply( $test_href, $expected_href, q{Leave hash unchanged} );

done_testing();
