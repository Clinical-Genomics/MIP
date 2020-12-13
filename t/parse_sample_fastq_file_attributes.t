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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::File_info}      => [qw{ parse_sample_fastq_file_attributes }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ parse_sample_fastq_file_attributes };

diag(   q{Test parse_sample_fastq_file_attributes from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $DATE        => q{161011};
Readonly my $READ_LENGTH => 151;

my $log = test_log( {} );

## Given fastq files
my %file_info = (
    ADM1059A1 => {
        mip_infiles_dir =>
          catdir( $Bin, qw{ data 643594-miptest test_data ADM1059A1 fastq } ),
        mip_infiles => [qw{ 1_161011_TestFilev2_ADM1059A1_TCCGGAGA_1.fastq.gz }],
    },
);
my $sample_id = q{ADM1059A1};

my %infile_info = parse_sample_fastq_file_attributes(
    {
        file_info_href => \%file_info,
        file_name      => $file_info{$sample_id}{mip_infiles}[0],
        infiles_dir    => $file_info{$sample_id}{mip_infiles_dir},
        sample_id      => $sample_id,
    }
);
my %expected_infile_info = (
    date              => $DATE,
    direction         => 1,
    flowcell          => q{TestFilev2},
    index             => q{TCCGGAGA},
    infile_sample_id  => q{ADM1059A1},
    is_interleaved    => undef,
    lane              => 1,
    read_file_command => q{gzip -d -c},
    read_length       => $READ_LENGTH,
);
my @expected_lanes = (1);

## Then file info for sample fastq should be collected
is_deeply( \%infile_info, \%expected_infile_info,
    q{Set and returned file info for sample fastq file} );

is_deeply( $file_info{$sample_id}{lanes},
    \@expected_lanes, q{Set sample lanes in file_info} );

done_testing();
