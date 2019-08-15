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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Program::Variantcalling::Bcftools} => [qw{ bcftools_rename_vcf_samples }],
        q{MIP::Test::Fixtures}                    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Bcftools qw{ bcftools_rename_vcf_samples };

diag(   q{Test bcftools_rename_vcf_samples from Bcftools.pm v}
      . $MIP::Program::Variantcalling::Bcftools::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given
my @sample_ids          = qw{ sample_1 };
my $test_dir            = File::Temp->newdir();
my $infile_path         = catfile(qw{ a infile.vcf });
my $outfile_path_prefix = q{outfile};

bcftools_rename_vcf_samples(
    {
        FILEHANDLE          => $FILEHANDLE,
        index               => 1,
        index_type          => q{csi},
        infile              => $infile_path,
        outfile_path_prefix => $outfile_path_prefix,
        output_type         => q{z},
        temp_directory      => catfile($test_dir),
        sample_ids_ref      => \@sample_ids,
    }
);

## Close the filehandle
close $FILEHANDLE;

## Then samle name file instruction should be written in file
my $create_sample_file = catfile( $test_dir, q{sample_name.txt} );
my ($is_ok_create_sample_file) = $file_content =~ /($create_sample_file)/xms;
ok( $is_ok_create_sample_file, q{Wrote reheader sample names file} );

## Then reheader instruction should be written in file
my ($is_ok_reheader) = $file_content =~ /(bcftools \s+ reheader)/xms;
ok( $is_ok_reheader, q{Wrote reheader instructions} );

## Then view instruction should be written in file
my ($is_ok_view) = $file_content =~ /(bcftools \s+ view)/xms;
ok( $is_ok_view, q{Wrote view and index instructions} );

done_testing();
