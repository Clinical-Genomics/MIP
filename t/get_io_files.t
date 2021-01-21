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

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::File_info} => [qw{ set_io_files }],
        q{MIP::Get::File} => [qw{ get_io_files }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_io_files };
use MIP::Get::File qw{ get_io_files };

diag(   q{Test get_io_files from File.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my $chain_main     = q{MAIN};
my $chain_chanjo   = q{CHSEX};
my $id             = q{sample_1};
my $id_2           = q{sample_2};
my $temp_directory = catfile(qw{ a temp dir });

my @base_paths = (
    catfile(qw{ a dir fastq_sample_1_1.fastq.gz }),
    catfile(qw{ a dir fastq_sample_1_2.fastq.gz }),
);
my @base_temp_paths = (
    catfile( $temp_directory, q{fastq_sample_1_1.fastq.gz} ),
    catfile( $temp_directory, q{fastq_sample_1_2.fastq.gz} ),
);
my @file_paths =
  ( catfile(qw{ a test bwa_mem file_1.txt}), catfile(qw{ a test bwa_mem file_2.txt}), );
my @file_paths_2 = (
    catfile(qw{ a test picard_mergesamfiles file_1.txt}),
    catfile(qw{ a test picard_mergesamfiles file_2.txt}),
);
my @delly_file_paths =
  ( catfile(qw{ a test delly_call dc_file_1.bcf}), catfile(qw{ a test delly_call dc_file_2.bcf}), );

my %file_info = (
    $id => {
        mip_infiles_dir      => catfile(qw{ a dir }),
        mip_infiles          => [qw{ fastq_sample_1_1.fastq.gz fastq_sample_1_2.fastq.gz }],
        base                 => { file_paths => \@base_paths, },
        base_temp            => { file_paths => \@base_temp_paths, },
        bwa_mem              => { file_paths => \@file_paths, },
        picard_mergesamfiles => { file_paths => \@file_paths_2, },
        delly_call           => { file_paths => \@delly_file_paths, },
    },
);
my @order_recipes =
  qw{ bwa_mem picard_mergesamfiles markduplicates gatk_baserecalibration chanjo_sexcheck cnvnator_ar delly_call delly_reformat sv_combinevariantcallsets };

my %parameter = (
    bwa_mem                   => { chain             => $chain_main, },
    chanjo_sexcheck           => { chain             => $chain_chanjo, },
    cnvnator_ar               => { chain             => q{CNVNATOR}, },
    delly_call                => { chain             => q{DELLY_CALL}, },
    delly_reformat            => { chain             => q{DELLY_CALL}, },
    cache                     => { order_recipes_ref => \@order_recipes, },
    markduplicates            => { chain             => $chain_main, },
    picard_mergesamfiles      => { chain             => $chain_main, },
    sv_combinevariantcallsets => { chain             => q{CHAIN_SV}, },
);

my $recipe_name = q{bwa_mem};
my $stream      = q{in};
my $stream_out  = q{out};

## Given no set infiles - inherit from base i.e MIP infiles
my %io = get_io_files(
    {
        id             => $id,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
        stream         => $stream,
        temp_directory => $temp_directory,
    }
);

## Then infile for new recipe should be returned
is_deeply(
    \@{ $io{$stream}{file_paths} },
    \@{ $file_info{$id}{base}{file_paths} },
    q{Got fastq infile features for sample_1}
);

is_deeply(
    \@{ $io{temp}{file_paths} },
    \@{ $file_info{$id}{base_temp}{file_paths} },
    q{Got fastq tempfile features for sample_1}
);

## Given of bwa_mem for $id and $id_2 for chain main - set io out
set_io_files(
    {
        chain_id       => $chain_main,
        id             => $id,
        file_paths_ref => \@file_paths,
        file_info_href => \%file_info,
        recipe_name    => $recipe_name,
        stream         => $stream_out,
        temp_directory => $temp_directory,
    }
);

# Set for $id_2
set_io_files(
    {
        chain_id       => $chain_main,
        id             => $id_2,
        file_paths_ref => \@file_paths,
        file_info_href => \%file_info,
        recipe_name    => $recipe_name,
        stream         => $stream_out,
    }
);

## Given new recipe
my $merge_sam = q{picard_mergesamfiles};

%io = get_io_files(
    {
        id             => $id,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
        recipe_name    => $merge_sam,
        stream         => $stream,
    }
);

## Then outfile for bwa_mem should be returned as infiles for merge_sam for sample_1
is_deeply(
    \@{ $file_info{$id}{bwa_mem}{file_paths} },
    \@{ $io{$stream}{file_paths} },
    q{Got bwa_mem infile features for merge_sam for sample_1}
);

## Given another id
my %io_id_2 = get_io_files(
    {
        id             => $id_2,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
        recipe_name    => $merge_sam,
        stream         => $stream,
    }
);

## Then outfile for bwa_mem should be returned as infiles for merge_sam sample_2
is_deeply(
    \@{ $file_info{$id}{bwa_mem}{file_paths} },
    \@{ $io_id_2{$stream}{file_paths} },
    q{Got bwa_mem infile features for merge_sam for sample_2}
);

## Set new outfiles for second recipe
set_io_files(
    {
        chain_id       => $chain_main,
        id             => $id,
        file_paths_ref => \@file_paths_2,
        file_info_href => \%file_info,
        recipe_name    => $merge_sam,
        stream         => $stream_out,
    }
);

%io = get_io_files(
    {
        id             => $id,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
        recipe_name    => $merge_sam,
        stream         => $stream_out,
    }
);

## Then outfile for merge_sam should be returned
is_deeply(
    \@{ $file_info{$id}{picard_mergesamfiles}{file_paths} },
    \@{ $io{$stream_out}{file_paths} },
    q{Got picard_mergesamfiles outfile features for merge_sam for sample_1}
);

## Given the first recipe in a chain that should inherit from MAIN
my $first_in_chain_recipe_name = q{chanjo_sexcheck};
%io = get_io_files(
    {
        id             => $id,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
        recipe_name    => $first_in_chain_recipe_name,
        stream         => $stream,
    }
);

is_deeply(
    \@{ $file_info{$id}{picard_mergesamfiles}{file_paths} },
    \@{ $io{$stream}{file_paths} },
    q{Got picard_mergesamfiles outfile features as infiles for chanjo_sexcheck for sample_1}
);

## Given a recipe downstream of PARALLEL chain and other chain
my $downstream_recipe = q{sv_combinevariantcallsets};

%io = get_io_files(
    {
        id             => $id,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
        recipe_name    => $downstream_recipe,
        stream         => $stream,
    }
);

is_deeply(
    \@{ $file_info{$id}{picard_mergesamfiles}{file_paths} },
    \@{ $io{$stream}{file_paths} },
q{Got picard_mergesamfiles outfile features as infiles for sv_combinevariantcallsets for sample_1}
);

## Given a recipe in a PARALLEL chain and other chain
my $delly_call = q{delly_call};

%io = get_io_files(
    {
        id             => $id,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
        recipe_name    => $delly_call,
        stream         => $stream,
    }
);

is_deeply(
    \@{ $file_info{$id}{picard_mergesamfiles}{file_paths} },
    \@{ $io{$stream}{file_paths} },
    q{Got picard_mergesamfiles outfile features as infiles for delly_call for sample_1}
);

## Set new outfiles for second recipe in PARALLEL self chain
set_io_files(
    {
        chain_id       => q{DELLY_CALL},
        id             => $id,
        file_paths_ref => \@delly_file_paths,
        file_info_href => \%file_info,
        recipe_name    => $delly_call,
        stream         => $stream_out,
    }
);

## Given a recipe in a PARALLEL chain and other chain
my $delly_reformat = q{delly_reformat};

%io = get_io_files(
    {
        id             => $id,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
        recipe_name    => $delly_reformat,
        stream         => $stream,
    }
);

is_deeply(
    \@{ $file_info{$id}{delly_call}{file_paths} },
    \@{ $io{$stream}{file_paths} },
    q{Got delly_call outfile features as infiles for delly_reformat for sample_1}
);

## Given a already processed recipe
my %io_bwa = get_io_files(
    {
        id             => $id,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
        recipe_name    => $recipe_name,
        stream         => $stream,
        temp_directory => $temp_directory,
    }
);

## Then infile for recipe should be returned
is_deeply(
    \@{ $io_bwa{$stream}{file_paths} },
    \@{ $file_info{$id}{base}{file_paths} },
    q{Got persistent infile features for recipe and sample_1}
);

done_testing();
