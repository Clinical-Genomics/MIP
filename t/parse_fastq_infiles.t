#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Fastq}          => [qw{ parse_fastq_infiles }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Fastq qw{ parse_fastq_infiles };

diag(   q{Test parse_fastq_infiles from Fastq.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given compressed file, when proper data
my %active_parameter = (
    analysis_type => {
        ADM1059A1 => q{wgs},
        ADM1059A2 => q{wgs},
    },
    sample_ids => [qw{ ADM1059A1 }],
);
my %file_info = (
    ADM1059A1 => {
        mip_infiles_dir =>
          catdir( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq } ),
        mip_infiles => [qw{ 1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz }],
    },
);
my %sample_info;

parse_fastq_infiles(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        sample_info_href      => \%sample_info,
    }
);

## Then return undef
is( $file_info{is_files_compressed}{ADM1059A1}, 1, q{All files compressed} );

## Given uncompressed file
push @{ $active_parameter{sample_ids} }, q{ADM1059A2};
push @{ $file_info{ADM1059A2}{mip_infiles} },
  qw{ 1_161011_TestFilev2_ADM1059A2_CGCTCATT_1.fastq ADM1059A2.fastq.gz };
$file_info{ADM1059A2}{mip_infiles_dir} =
  catdir( $Bin, qw{ data 643594-miptest test_data bad_input } );

parse_fastq_infiles(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        sample_info_href      => \%sample_info,
    }
);

## Then return true
is( $file_info{is_files_compressed}{ADM1059A2},
    0, q{Files uncompressed and got run info from headers} );

## Remove ADM1059A2 from processing
pop @{ $active_parameter{sample_ids} };

## Given file, when no sample_id in file name
push @{ $active_parameter{sample_ids} }, q{ADM1059A3};
push @{ $file_info{ADM1059A3}{mip_infiles} }, qw{ 1_161011_TestFilev2_CGCTCATT_1.fastq };
$file_info{ADM1059A3}{mip_infiles_dir} =
  catdir( $Bin, qw{ data 643594-miptest test_data bad_input } );

trap {
    parse_fastq_infiles(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            sample_info_href      => \%sample_info,
        }
    );
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if sample_id in file name cannot be found} );
like(
    $trap->stderr,
    qr/Please \s+ check \s+ that \s+ the \s+ file \s+ name/xms,
    q{Throw fatal log message if sample_id in file name cannot be found}
);

done_testing();
