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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA        => q{,};
Readonly my $GRCH_VERSION => 37;
Readonly my $HG_VERSION   => 19;
Readonly my $SPACE        => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::File_info}      => [qw{ set_human_genome_reference_features }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_human_genome_reference_features };
use MIP::Log::MIP_log4perl qw{ initiate_logger };

diag(   q{Test set_human_genome_reference_features from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

my %parameter;

## Test hash
my %file_info;

## Test Ensemble genome
my $human_genome_reference = q{grch37_homo_sapiens_-d5-.fasta};
set_human_genome_reference_features(
    {
        file_info_href         => \%file_info,
        human_genome_reference => $human_genome_reference,
        parameter_href         => \%parameter,
    }
);

is( $file_info{human_genome_reference_version}, $GRCH_VERSION, q{grch version test} );
is( $file_info{human_genome_reference_source},  q{grch},       q{grch source test} );
is( $file_info{human_genome_reference_name_prefix},
    q{grch37_homo_sapiens_-d5-}, q{grch prefix test} );
is( $file_info{human_genome_compressed}, 0, q{grch compressed test} );
is( $parameter{human_genome_reference_file_endings}{build_file},
    undef, q{Did not set build file for compression } );

## Test Refseq genome
$human_genome_reference = q{hg19_homo_sapiens.fasta.gz};
set_human_genome_reference_features(
    {
        file_info_href         => \%file_info,
        human_genome_reference => $human_genome_reference,
        parameter_href         => \%parameter,
    }
);
is( $file_info{human_genome_reference_version}, $HG_VERSION, q{hg version test} );
is( $file_info{human_genome_reference_source},  q{hg},       q{hg source test} );
is( $file_info{human_genome_reference_name_prefix},
    q{hg19_homo_sapiens}, q{hg prefix test} );
is( $file_info{human_genome_compressed}, 1, q{hg compressed test} );
is( $parameter{human_genome_reference_file_endings}{build_file},
    1, q{Set build file for compression } );

## Given a human reference, when no version
# Clear for new test
%file_info = ();
my $human_genome_reference_no_version = q{grch_homo_sapiens_-d5-.fasta};
trap {
    set_human_genome_reference_features(
        {
            file_info_href         => \%file_info,
            human_genome_reference => $human_genome_reference_no_version,
            parameter_href         => \%parameter,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the version cannot be found} );
like(
    $trap->stderr,
    qr/MIP \s+ cannot \s+ detect/xms,
    q{Throw fatal log message if the version cannot be found}
);

done_testing();
