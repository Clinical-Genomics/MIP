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
use MIP::Constants qw{ $COMMA $EMPTY_STR $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Download}       => [qw{ check_user_reference }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Download qw{ check_user_reference };

diag(   q{Test check_user_reference from Download.pm v}
      . $MIP::Download::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => uc q{mip_download}, } );

## Given some references
my %active_parameter =
  test_mip_hashes( { mip_hash_name => q{download_active_parameter}, } );

## Clean-up since there are already references in there
delete $active_parameter{reference};

## When using a bad reference key
$active_parameter{reference}{not_a_reference} = [qw{ decoy_5 }];

trap {
    check_user_reference(
        {
            reference_genome_versions_ref =>
              \@{ $active_parameter{reference_genome_versions} },
            reference_ref               => \%{ $active_parameter{reference_feature} },
            user_supplied_reference_ref => \%{ $active_parameter{reference} },
        }
    )
};

## Then exit and throw warn log message
ok( $trap->exit, q{Exit if the reference key cannot be found} );
like(
    $trap->stderr,
    qr/Cannot \s+ find \s+ reference \s+ key/xms,
    q{Throw warn log message if the reference key cannot be found}
);

# Clean-up
delete $active_parameter{reference};

## Given a valid reference key

## When using a bad reference version
$active_parameter{reference}{human_reference} = [qw{ bad_version }];

trap {
    check_user_reference(
        {
            reference_genome_versions_ref =>
              \@{ $active_parameter{reference_genome_versions} },
            reference_ref               => \%{ $active_parameter{reference_feature} },
            user_supplied_reference_ref => \%{ $active_parameter{reference} },
        }
    )
};

## Then throw warn log message
like(
    $trap->stderr,
    qr/Cannot \s+ find \s+ version \s+ key/xms,
    q{Throw warn log message if the reference version cannot be found}
);

## When using a valid reference
$active_parameter{reference}{human_reference} = [qw{ decoy_5 }];

my @returns = trap {
    check_user_reference(
        {
            reference_genome_versions_ref =>
              \@{ $active_parameter{reference_genome_versions} },
            reference_ref               => \%{ $active_parameter{reference_feature} },
            user_supplied_reference_ref => \%{ $active_parameter{reference} },
        }
    )
};
## Then return true
ok( $returns[0], q{Checked reference} );

## Then throw no warning
is( $trap->stderr, $EMPTY_STR, q{Found all references} );

done_testing();
