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
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::File::Decompression} => [qw{ decompress_file }],
        q{MIP::Test::Fixtures}      => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Decompression qw{ decompress_file };

diag(   q{Test decompress_file from Decompression.pm v}
      . $MIP::File::Decompression::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

## Given an already decompressed option
my $decompress_program;
my $outdir_path  = catdir(qw{ a dir });
my $outfile_path = catfile(qw{ a dir file.fastq.gz});

my $return = decompress_file(
    {
        FILEHANDLE         => $FILEHANDLE,
        outdir_path        => $outdir_path,
        outfile_path       => $outfile_path,
        decompress_program => $decompress_program,
    }
);

## Then return undef
is( $return, undef, q{Already decompressed file} );

# For storing info to write
my $file_content;

my @programs = qw{ gzip unzip tar };

PROGRAM:
foreach my $program (@programs) {

## Store file content in memory by using referenced variable
    open $FILEHANDLE, q{>}, \$file_content
      or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given a compressed file
    decompress_file(
        {
            FILEHANDLE         => $FILEHANDLE,
            outdir_path        => $outdir_path,
            outfile_path       => $outfile_path,
            decompress_program => $program,
        }
    );

## Close the filehandle
    close $FILEHANDLE;

    my ($returned_base_command) = $file_content =~ /^($program)/xms;

## Then the base command should be in file content
    is( $returned_base_command, $program, q{ Decompressed using } . $program );
}

done_testing();
