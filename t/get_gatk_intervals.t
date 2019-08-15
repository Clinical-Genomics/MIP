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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COLON => q{:};
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Get::Parameter} => [qw{ get_gatk_intervals }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_gatk_intervals };

diag(   q{Test get_gatk_intervals from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

my %expected_output = (
    wes => {
        1 => [ catfile(qw{ a dir 1_bed_file.bed.interval_list }) ],
        2 => [ catfile(qw{ a dir 2_bed_file.bed.interval_list }) ],
    },
    wgs => {
        1 => [1],
        2 => [2],
    },
    wts => {
        1 => [1],
        2 => [2],
    },
);

my %exome_target_bed = ( q{bed_file.bed} => q{test_sample}, );

ANALYSIS_TYPE:
foreach my $analysis_type ( keys %expected_output ) {

    ## Given the input below
    my %gatk_intervals = get_gatk_intervals(
        {
            analysis_type         => $analysis_type,
            contigs_ref           => [qw{ 1 2 }],
            FILEHANDLE            => $FILEHANDLE,
            outdirectory          => catdir(qw{ a dir }),
            reference_dir         => catdir(qw{ a dir reference_dir}),
            exome_target_bed_href => \%exome_target_bed,
            file_ending           => q{.interval_list},
            log                   => $log,
            sample_id             => q{test_sample},
        }
    );

    ## Then, generate gatk intervals. Chromosomes for WGS/WTS and paths to contig_bed_files for WES
    is_deeply(
        \%gatk_intervals,
        $expected_output{$analysis_type},
        qq{Set gatk intervals for $analysis_type}
    );
}

close $FILEHANDLE;

done_testing();
