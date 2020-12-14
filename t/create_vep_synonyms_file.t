#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
        q{MIP::File::Format::Vep} => [qw{ create_vep_synonyms_file }],
        q{MIP::Test::Fixtures}    => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Vep qw{ create_vep_synonyms_file };

diag(   q{Test create_vep_synonyms_file from Vep.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log          = test_log( { no_screen => 1, } );
my $test_dir     = File::Temp->newdir();
my $outfile_path = catfile( $test_dir, q{synonyms.tsv} );

## Given
my $bad_version = q{not_valid_version};

my $is_ok = create_vep_synonyms_file(
    {
        log          => $log,
        outfile_path => $outfile_path,
        version      => $bad_version,
    }
);

## Then return false
is( $is_ok, undef, q{Return if no defined synonyms map} );

## Given
my $valid_version = q{38};

my $ret_outfile_path = create_vep_synonyms_file(
    {
        log          => $log,
        outfile_path => $outfile_path,
        version      => $valid_version,
    }
);

## Then return outfile_path for synonyms file
is( $outfile_path, $ret_outfile_path, q{Return outfile_path if defined synonyms map} );
ok( -e $outfile_path, q{Created synonyms file} );
done_testing();
