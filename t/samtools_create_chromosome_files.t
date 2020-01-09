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
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Program::Samtools} => [qw{ samtools_create_chromosome_files }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Samtools qw{ samtools_create_chromosome_files };

diag(   q{Test samtools_create_chromosome_files from Samtools.pm v}
      . $MIP::Program::Samtools::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $filehandle = IO::Handle->new();

# For storing info to write
my $file_content;

# Store file content in memory by using referenced variable
open $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given request to index
samtools_create_chromosome_files(
    {
        filehandle         => $filehandle,
        infile_path        => catfile(qw{ path test.fa}),
        max_process_number => 1,
        regions_ref        => [qw{ X Y }],
        suffix             => q{.fai},
        temp_directory     => catdir(qw{ temp }),
    }
);

# Close the filehandle
close $filehandle;

## Then split and index each region
my @faidx_regexps = (
    qr{samtools\sfaidx\spath\/test\.fa\sX\s>\stemp\/X\.fai}xms,
    qr{wait}xms, qr{samtools\sfaidx\spath\/test\.fa\sY\s>\stemp\/Y\.fai}xms,
);

my $iterator;
REGEX:
foreach my $regex (@faidx_regexps) {
    $iterator++;
    my ($faidx_command) = $file_content =~ $regex;
    ok( $faidx_command, qq{Test command $iterator} );
}

done_testing();
