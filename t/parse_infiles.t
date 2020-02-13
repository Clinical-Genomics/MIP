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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Parse::Parameter} => [qw{ parse_infiles }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Parameter qw{ parse_infiles };

diag(   q{Test parse_infiles from Parameter.pm v}
      . $MIP::Parse::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given proper test data, when parsing
my $sample_id        = q{ADM1059A1};
my %active_parameter = (
    infile_dirs => {
        catfile( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq } ) =>
          q{ADM1059A1},
    },
    sample_ids => [qw{ ADM1059A1 }],
);
my %file_info;

my @returns = trap {
    parse_infiles(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            log                   => $log,
        }
    )
};

## Then return true
is( $returns[0], 1, q{Parsed input files} );

## Then broacast info on infile and indir per sample_id
like( $trap->stderr, qr/INFO/xms, q{Broadcast INFO log message of parsing} );

## Given no infile_directory
delete $active_parameter{infile_dirs};

trap {
    parse_infiles(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            log                   => $log,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the infile_dirs parameter cannot be found} );
like( $trap->stderr, qr/Could\s+not\s+detect/xms,
    q{Throw fatal log message if the infile_dirs parameter cannot be found} );

done_testing();
