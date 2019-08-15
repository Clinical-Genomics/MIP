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

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given a file to decompress, when method gzip unzip and tar
my $outfile_path  = q{a_file_path.gz};
my $reference_dir = q{a_reference_dir};
my @programs      = qw{ gzip unzip tar };

PROGRAM:
foreach my $program (@programs) {

    decompress_file(
        {
            FILEHANDLE         => $FILEHANDLE,
            outdir_path        => $reference_dir,
            outfile_path       => $outfile_path,
            decompress_program => $program,
        }
    );

}
## Close the filehandle
close $FILEHANDLE;

PROGRAM:
foreach my $program (@programs) {

## Then each command should be in the same file
    my ($returned_base_command) = $file_content =~ /^($program)/xms;
    ok( $returned_base_command, qq{Wrote $program command for decompression} );

}

done_testing();
