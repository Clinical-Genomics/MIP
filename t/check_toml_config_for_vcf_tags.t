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
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_constants test_log };
use MIP::Test::Writefile qw{ write_toml_config };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Reference}      => [qw{ check_toml_config_for_vcf_tags }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ check_toml_config_for_vcf_tags };

diag(   q{Test check_toml_config_for_vcf_tags from Reference.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create log object
test_log( { log_name => q{MIP}, no_screen => 0, } );

my $cluster_reference_path = catdir( dirname($Bin), qw{ t data references } );
my $toml_template_path =
  catfile( $cluster_reference_path, q{grch37_vcfanno_config_template-v1.0-.toml} );
my $toml_config_path =
  catfile( $cluster_reference_path, q{grch37_vcfanno_config-v1.0-.toml} );

## Update path in toml config
write_toml_config(
    {
        test_reference_path => $cluster_reference_path,
        toml_config_path    => $toml_config_path,
        toml_template_path  => $toml_template_path,
    }
);
## Given a toml config with all vcf tags present in vcf
my %active_parameter_test = (
    variant_annotation => 1,
    vcfanno_config     => $toml_config_path,
);
my %test_process_return = (
    buffers_ref   => [],
    error_message => undef,
    stderrs_ref   => [],
    stdouts_ref   => [q{AF AF_POPMAX}],
    success       => 1,
);
test_constants(
    {
        test_process_return_href => \%test_process_return,
    }
);

## Then all is ok
my $is_ok = check_toml_config_for_vcf_tags(
    {
        active_parameter_href => \%active_parameter_test,
    }
);
ok( $is_ok, q{Vcf file contains required annotations} );

## Given a vcfanno annotation request that tries to use a non-existing vcf tag
$toml_template_path =
  catfile( $cluster_reference_path, q{grch37_vcfanno_config_bad_template-v1.0-.toml} );
$toml_config_path =
  catfile( $cluster_reference_path, q{grch37_vcfanno_config-v1.0-.toml} );

write_toml_config(
    {
        test_reference_path => $cluster_reference_path,
        toml_config_path    => $toml_config_path,
        toml_template_path  => $toml_template_path,
    }
);
$test_process_return{stdouts_ref} = [];
test_constants(
    {
        test_process_return_href => \%test_process_return,
    }
);

trap {
    check_toml_config_for_vcf_tags(
        {
            active_parameter_href => \%active_parameter_test,
        }
    )
};

## Then print fatal log message and exit
ok( $trap->exit, q{Exit when vcf annotation tags are missing } );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

unlink $toml_config_path;
done_testing();
