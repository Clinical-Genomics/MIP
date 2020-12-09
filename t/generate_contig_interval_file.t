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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Contigs}        => [qw{ generate_contig_interval_file }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Contigs qw{ generate_contig_interval_file };

diag(   q{Test generate_contig_interval_file from Contigs.pm v}
      . $MIP::Contigs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $file_content;

## Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given
my @contigs               = qw{ 1 2 };
my $exome_target_bed_file = q{grch37_agilent_sureselect_targets_-v5-.bed};
my $reference_dir         = catfile( $Bin, qw{ data reference } );
my $outdirectory          = catfile(qw{ a outdir });

my %bed_file_path = generate_contig_interval_file(
    {
        contigs_ref           => \@contigs,
        exome_target_bed_file => $exome_target_bed_file,
        filehandle            => $filehandle,
        max_cores_per_node    => 1,
        outdirectory          => $outdirectory,
        reference_dir         => $reference_dir,
    }
);
my %expected_bed_file_path = (
    1 => [ catfile( $outdirectory, q{1_} . $exome_target_bed_file ) ],
    2 => [ catfile( $outdirectory, q{2_} . $exome_target_bed_file ) ]
);

## Then contig bed file paths for each contig should be returned
is_deeply( \%bed_file_path, \%expected_bed_file_path, q{Generated bed file} );

## Given a file ending
my $file_suffix = q{.custom};

## When writing instructions to generate seperate contig files
my %bed_file_path_with_ending = generate_contig_interval_file(
    {
        contigs_ref           => \@contigs,
        exome_target_bed_file => $exome_target_bed_file,
        filehandle            => $filehandle,
        file_ending           => $file_suffix,
        max_cores_per_node    => 1,
        outdirectory          => $outdirectory,
        reference_dir         => $reference_dir,
    }
);
my %expected_bed_file_path_with_ending = (
    1 => [ catfile( $outdirectory, q{1_} . $exome_target_bed_file . $file_suffix ) ],
    2 => [ catfile( $outdirectory, q{2_} . $exome_target_bed_file . $file_suffix ) ]
);
## Close the filehandle
close $filehandle;

## Then contig bed file paths with file suffix for each contig should be returned
is_deeply(
    \%bed_file_path_with_ending,
    \%expected_bed_file_path_with_ending,
    q{Generated bed file with supplied file suffix}
);

done_testing();
