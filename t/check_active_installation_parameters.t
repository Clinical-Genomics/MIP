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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Mip_install}    => [qw{ check_active_installation_parameters }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Mip_install qw{ check_active_installation_parameters };

diag(   q{Test check_active_installation_parameters from Mip_install.pm v}
      . $MIP::Mip_install::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given sbatch mode and a project id
my $project_id  = q{test};
my $sbatch_mode = 1;

## When installing in sbatch mode with a project id
my $is_ok = check_active_installation_parameters(
    {
        project_id  => $project_id,
        sbatch_mode => $sbatch_mode,
    }
);

## Then return true
ok( $is_ok, q{Project id was supplied in sbatch mode} );

## Given sbatch mode

## When installing in sbatch mode with no project id
trap {
    check_active_installation_parameters(
        {
            project_id  => undef,
            sbatch_mode => $sbatch_mode,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if no project id in sbatch mode} );
like(
    $trap->stderr,
    qr/The \s+ option \s+  "project_id"/xms,
    q{Throw fatal log message if no project id in sbatch mode}
);
done_testing();
