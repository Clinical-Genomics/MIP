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
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_log test_standard_cli };

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
        q{MIP::Check::Download} => [qw{ check_user_reference }],
        q{MIP::Test::Fixtures}  => [qw{ test_mip_hashes test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Download qw{ check_user_reference };

diag(   q{Test check_user_reference from Download.pm v}
      . $MIP::Check::Download::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { log_name => uc q{mip_download}, } );

## Given a bad reference key
my %active_parameter =
  test_mip_hashes( { mip_hash_name => q{download_active_parameter}, } );

## Clean-up since there are already references in there
delete $active_parameter{reference};
$active_parameter{reference}{not_a_reference} = [qw{ decoy_5 }];

trap {
    check_user_reference(
        {
            user_supplied_reference_ref => \%{ $active_parameter{reference} },
            reference_genome_versions_ref =>
              \@{ $active_parameter{reference_genome_versions} },
            reference_ref => \%{ $active_parameter{reference_feature} },
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

## Given a valid reference key with a bad reference version

$active_parameter{reference}{human_reference} = [qw{ bad_version }];

trap {
    check_user_reference(
        {
            user_supplied_reference_ref => \%{ $active_parameter{reference} },
            reference_genome_versions_ref =>
              \@{ $active_parameter{reference_genome_versions} },
            reference_ref => \%{ $active_parameter{reference_feature} },
        }
    )
};

## Then throw warn log message
like(
    $trap->stderr,
    qr/Cannot \s+ find \s+ version \s+ key/xms,
    q{Throw warn log message if the reference key cannot be found}
);

## Given an valid reference
$active_parameter{reference}{human_reference} = [qw{ decoy_5 }];

my $is_ok = check_user_reference(
    {
        user_supplied_reference_ref => \%{ $active_parameter{reference} },
        reference_genome_versions_ref =>
          \@{ $active_parameter{reference_genome_versions} },
        reference_ref => \%{ $active_parameter{reference_feature} },
    }
);

## Then return true
ok( $is_ok, q{Checked reference} );

done_testing();
