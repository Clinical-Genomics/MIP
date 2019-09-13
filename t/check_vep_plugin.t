#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Check::Parameter} => [qw{ check_vep_plugin }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Parameter qw{ check_vep_plugin };

diag(   q{Test check_vep_plugin from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given an undefine vep plugin hash
my %vep_plugin;

my $is_not_ok = check_vep_plugin(
    {
        log             => $log,
        parameter_name  => q{vep_plugin},
        vep_plugin_href => \%vep_plugin,
    }
);
## Then there is nothing to check - return false
is( $is_not_ok, 0, q{Nothing to check} );

%vep_plugin = (
    dbNSFP => {
        exists_check => q{directory},
        path         => cwd(),
        parameters   => [qw{ param_1 param_2 }],
    },
);

my $is_ok = check_vep_plugin(
    {
        log             => $log,
        parameter_name  => q{vep_plugin},
        vep_plugin_href => \%vep_plugin,
    }
);

## Then
ok( $is_ok, q{Checked vep plugin hash} );

## Given a not valid hash ref
$vep_plugin{not_valid_annotation} = [q{not_a_valid_ref}];

trap {
    check_vep_plugin(
        {
            log             => $log,
            parameter_name  => q{vep_plugin},
            vep_plugin_href => \%vep_plugin,
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if not a hash ref } );
like( $trap->die, qr/Is\s+not\s+a/xms, q{Not a hash ref} );

delete $vep_plugin{not_valid_annotation};

done_testing();
