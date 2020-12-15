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
        q{MIP::Reference}      => [qw{ write_contigs_size_file }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ write_contigs_size_file };

diag(   q{Test write_contigs_size_file from Reference.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

my $temp_dir = File::Temp->newdir();

## Given proper data
my $fai_file_path =
  catfile( $Bin, qw{ data references grch37_homo_sapiens_-d5-.fasta.fai } );
my $outfile_path = catfile( $temp_dir, q{chromosome_size_file.tsv} );

write_contigs_size_file(
    {
        fai_file_path => $fai_file_path,
        outfile_path  => $outfile_path,
    }
);

## Then write the fai contig names and lengths
ok( -e $outfile_path, q{Wrote contig file size} );

## Given a file that does not exist
my $wrong_file = catfile( $Bin, q{chromosome_size_does_not_exist} );

trap {
    write_contigs_size_file(
        {
            fai_file_path => catfile( $Bin, q{chromosome_size_does_not_exist} ),
            outfile_path  => $outfile_path,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if contigs cannot be written} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if contigs cannot be found} );

done_testing();
