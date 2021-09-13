#! /usr/bin/env perl

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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Environment::Container} => [qw{ check_installed_containers }],
        q{MIP::Test::Fixtures}         => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Container qw{ check_installed_containers };

diag(   q{Test check_installed_containers from Container.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { no_screen => 0, } );

## Given docker hub address and executable container
my %container = (
    fastqc => {
        executable => {
            fastqc => undef,
        },
        uri => catfile(
            dirname($Bin), qw{ t data modules miniconda envs mip_ci bin fastqc:0.11.9.sif }
        ),
    },
    htslib => {
        executable => {
            bgzip    => undef,
            samtools => undef,
            tabix    => undef,
        },
        uri => q{docker.io/clinicalgenomics/htslib:1.13},
    },
);

## Then everything is OK
my $ok = check_installed_containers(
    {
        container_href => \%container,
    }
);
ok( $ok, q{Sif's exists and are executable} );

## Given missing executable
%container = (
    manta => {
        executable => {
            q{configManta.py} => undef,
            q{runWorkflow.py} => q{no_executable_in_image},
        },
        uri => catfile(qw{ missing file manta:1.6.0.sif }),
    },
);

## Then Croak with error message
trap {
    check_installed_containers(
        {
            container_href => \%container,
        }
    )
};

ok( $trap->exit, q{Exit if sif cannot be found} );
like(
    $trap->stderr,
    qr/Could\snot\sfind\sintended\smanta\sexecutable_file/xms,
    q{Print log message for missing sif}
);

## Given non executable sif
%container = (
    multiqc => {
        executable => {
            multiqc => undef,
        },
        uri =>
          catfile( dirname($Bin), qw{t data modules miniconda envs mip_ci bin multiqc:v1.11.sif } ),
    },
);

## Then Croak with error message
trap {
    check_installed_containers(
        {
            container_href => \%container,
        }
    )
};

ok( $trap->exit, q{Exit if sif isn't executable} );
like(
    $trap->stderr,
    qr/Could\snot\sfind\sintended\smultiqc\sexecutable_file/xms,
    q{Print log message for non executable sif}
);

done_testing();
