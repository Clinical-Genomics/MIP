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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE $NEWLINE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::File::Decompression} => [qw{ decompress_files }],
        q{MIP::Test::Fixtures}      => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Decompression qw{ decompress_files };

diag(   q{Test decompress_files from Decompression.pm v}
      . $MIP::File::Decompression::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a file to decompress, when method gzip unzip and tar
my $file_paths_ref = [q{a_file_path.gz}];
my $file_path      = q{a_file_path.gz};
my $outdir_path    = q{a_outdir};

my %program = (
    gzip  => q{a_outfile.gz},
    tar   => q{a_outfile.tar},
    unzip => q{a_outfile.zip},
);

PROGRAM:
while ( my ( $program, $outfile_path ) = each %program ) {

    my @decompress_commands = decompress_files(
        {
            file_path    => $file_path,
            outdir_path  => $outdir_path,
            outfile_path => $outfile_path,
            program      => $program,
        }
    );
    ## Then write matching command
    is( $decompress_commands[0], $program, qq{Wrote $program command for decompression} );
}

## Given missing file_paths for tar
my $outfile_path = q{a_outfile};

trap {
    decompress_files(
        {
            file_paths_ref => $file_paths_ref,
            outdir_path    => $outdir_path,
            outfile_path   => $outfile_path,
            program        => q{tar},
        }
    )
};
## Then croak
ok( $trap->die, q{Croak when missing infile} );

done_testing();
