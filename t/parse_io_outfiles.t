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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::File_info} => [qw{ parse_io_outfiles }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ parse_io_outfiles };

diag(   q{Test parse_io_outfiles from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a chain id
my $chain_main = q{CHAIN_MAIN};

## Given an id
my $id = q{sample_1};

## Given base file paths for an id in file_info hash
my %file_info = (
    $id => {
        mip_infiles_dir => catfile(qw{ a dir }),
        mip_infiles     => [qw{ fastq_sample_1_1.fastq.gz fastq_sample_1_2.fastq.gz }],
    },
);

## Given an order of recipes
my @order_recipes =
  qw{ bwa_mem picard_mergesamfiles markduplicates gatk_baserecalibration chanjo_sexcheck cnvnator_ar sv_combinevariantcallsets };

## Given a parameter hash with recipe names, chains and a recipe order
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

## Given a recipe name
my $recipe_name = q{bwa_mem};

## Given no set infiles

## When parsing io outfile
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

## Then inherit outfiles from base i.e MIP infiles for recipe
is_deeply(
    \%{ $file_info{io}{$chain_main}{$id}{$recipe_name}{out} },
    \%{ $io{out} },
    q{Set and got fastq file features for sample_1}
);

## Given set of iterators and infile prefix - construct default paths with iterator

## When parsing io outfile
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
