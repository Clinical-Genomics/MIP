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
        q{MIP::Validate::Case} => [qw{ check_infiles }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Validate::Case qw{ check_infiles };

diag(   q{Test check_infiles from Case.pm v}
      . $MIP::Validate::Case::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Proper data
my @infiles          = qw{file_1_sample-1.fastq.gz file_2_sample-1.fastq.gz};
my $infile_directory = q{a_test_dir};
my $sample_id        = q{sample-1};

## Given no infile
my @no_infiles;

trap {
    check_infiles(
        {
            infiles_ref      => \@no_infiles,
            infile_directory => $infile_directory,
            sample_id        => $sample_id,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if no infiles} );
like(
    $trap->stderr,
    qr/Could \s+ not \s+ find \s+ any \s+ fastq/xms,
    q{Throw fatal log message if no infiles}
);

## Given incorrect sample id
my $wrong_sample_id = q{wrong_sample_id};

trap {
    check_infiles(
        {
            infiles_ref      => \@infiles,
            infile_directory => $infile_directory,
            sample_id        => $wrong_sample_id,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if wrong sample id} );
like(
    $trap->stderr,
    qr/Could \s+ not \s+ detect \s+ sample_id/xms,
    q{Throw fatal log message if wrong sample id}
);

## Given wrong sample_id in infile name
my @wrong_infiles = qw{file_1_wrong_sample_id.fastq.gz file_2_sample-1.fastq.gz};

trap {
    check_infiles(
        {
            infiles_ref      => \@wrong_infiles,
            infile_directory => $infile_directory,
            sample_id        => $sample_id,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if wrong sample id in file name} );
like(
    $trap->stderr,
    qr/Check \s+ that:/xms,
    q{Throw fatal log message if wrong sample id in file name}
);

## Given proper indata
my $is_ok = check_infiles(
    {
        infiles_ref      => \@infiles,
        infile_directory => $infile_directory,
        sample_id        => $sample_id,
    }
);

## Then return true
ok( $is_ok, q{Passed check} );

done_testing();
