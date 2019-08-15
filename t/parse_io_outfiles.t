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
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Parse::File}    => [qw{ parse_io_outfiles }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::File qw{ parse_io_outfiles };

diag(   q{Test parse_io_outfiles from File.pm v}
      . $MIP::Parse::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my $chain_main = q{CHAIN_MAIN};
my $id         = q{sample_1};
my @file_paths =
  ( catfile(qw{ a test bwa_mem file_1.txt}), catfile(qw{ a test bwa_mem file_2.txt}), );
my %file_info = (
    $id => {
        mip_infiles_dir => catfile(qw{ a dir }),
        mip_infiles     => [qw{ fastq_sample_1_1.fastq.gz fastq_sample_1_2.fastq.gz }],
    },
);
my @order_recipes =
  qw{ bwa_mem picard_mergesamfiles markduplicates gatk_baserecalibration chanjo_sexcheck cnvnator_ar sv_combinevariantcallsets };

my %parameter = (
    bwa_mem => {
        chain          => $chain_main,
        outfile_suffix => q{.bam},
    },
    chanjo_sexcheck           => { chain             => q{CHAIN_CHANJO}, },
    cnvnator_ar               => { chain             => q{CNVNATOR}, },
    cache                     => { order_recipes_ref => \@order_recipes, },
    markduplicates            => { chain             => $chain_main, },
    picard_mergesamfiles      => { chain             => $chain_main, },
    sv_combinevariantcallsets => { chain             => q{CHAIN_SV}, },
);
my $recipe_name = q{bwa_mem};

## Given no set infiles - inherit from base i.e MIP infiles
my %io = parse_io_outfiles(
    {
        chain_id       => $chain_main,
        id             => $id,
        file_info_href => \%file_info,
        file_paths_ref => [ catfile(qw{ a dir file.bam}) ],
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
    }
);

## Then outfile for recipe should be returned
is_deeply(
    \%{ $file_info{io}{$chain_main}{$id}{$recipe_name}{out} },
    \%{ $io{out} },
    q{Set and got fastq file features for sample_1}
);

## Given set of iterators and infile prefix - construct default paths with iterator
%io = parse_io_outfiles(
    {
        chain_id         => $chain_main,
        file_info_href   => \%file_info,
        file_name_prefix => q{alignment_file},
        id               => $id,
        iterators_ref    => [qw{ MT X }],
        outdata_dir      => catdir(qw{ test dir }),
        parameter_href   => \%parameter,
        recipe_name      => $recipe_name,
    }
);

## Then construct outfile paths
my %expected_outfile_hash = (
    MT => catfile( qw{test dir}, $id, $recipe_name, q{alignment_file.MT.bam} ),
    X  => catfile( qw{test dir}, $id, $recipe_name, q{alignment_file.X.bam} ),
);
is_deeply( \%{ $io{out}{file_path_href} },
    \%expected_outfile_hash, q{Construct default filenames from prefix and iterators} );

done_testing();
