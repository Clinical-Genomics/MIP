#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
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
        q{MIP::File_info}      => [qw{ set_infiles }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_infiles };

diag(   q{Test set_infiles from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given infile and indir to set
my $sample_id = q{sample-1};
my %file_info;
my $infile_directory = q{a_test_dir};
my @infiles          = qw{file_1_sample-1.fastq.gz file_2_samnple-1.fastq.gz};

set_infiles(
    {
        file_info_href   => \%file_info,
        infiles_ref      => \@infiles,
        infile_directory => $infile_directory,
        sample_id        => $sample_id,
    }
);
my %expected_indir_path = ( q{sample-1} => { mip_infiles_dir => $infile_directory, }, );
my %expected_infile     = ( q{sample-1} => { mip_infiles     => [@infiles], }, );

## Then set it in the infile and indir hash
is(
    $file_info{q{sample-1}}{mip_infiles_dir},
    $expected_indir_path{q{sample-1}}{mip_infiles_dir},
    q{Set indir_path hash}
);

is_deeply(
    \@{ $file_info{q{sample-1}}{mip_infiles} },
    \@{ $expected_infile{q{sample-1}}{mip_infiles} },
    q{Set infile hash}
);

done_testing();
