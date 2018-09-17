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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Get::File}      => [qw{ get_io_files }],
        q{MIP::Set::File}      => [qw{ set_io_files }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::File qw{ get_io_files };
use MIP::Set::File qw{ set_io_files };

diag(   q{Test get_io_files from File.pm v}
      . $MIP::Get::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my $chain_main   = q{CHAIN_MAIN};
my $chain_chanjo = q{CHAIN_CHSEX};
my $id           = q{sample_1};
my $id_2         = q{sample_2};
my @file_paths   = (
    catfile(qw{ a test bwa_mem file_1.txt}),
    catfile(qw{ a test bwa_mem file_2.txt}),
);
my %file_info = (
    $id => {
        mip_infiles => [
            catfile(qw{ a dir fastq_sample_1_1.fastq.gz }),
            catfile(qw{ a dir fastq_sample_1_2.fastq.gz })
        ],
    },
);
my @order_programs =
  qw{ bwa_mem picard_mergesamfiles markduplicates gatk_baserecalibration chanjo_sexcheck cnvnator};

my %parameter = (
    bwa_mem                   => { chain => $chain_main, },
    picard_mergesamfiles      => { chain => $chain_main, },
    markduplicates            => { chain => $chain_main, },
    chanjo_sexcheck           => { chain => $chain_chanjo, },
    cnvnator                  => { chain => q{CNVNATOR}, },
    sv_combinevariantcallsets => { chain => q{CHAIN_SV}, },
);
my $program_name = q{picard_mergesamfiles};

## Given no set files - inherit from base i.e MIP infiles
my %io = get_io_files(
    {
        chain_id           => $chain_main,
        id                 => $id,
        file_info_href     => \%file_info,
        order_programs_ref => \@order_programs,
        parameter_href     => \%parameter,
        program_name       => $program_name,
    }
);

## Then infile for new program should be returned
is_deeply( \%{ $file_info{io}{$chain_main}{$id} },
    \%io, q{Got fastq file features for sample_1} );

#use Data::Printer;
#p(%io);

## Set for $id
## Program name does not matter - features gets overwritten
set_io_files(
    {
        chain_id       => $chain_main,
        id             => $id,
        file_paths_ref => \@file_paths,
        file_info_href => \%file_info,
    }
);

## Set for $id_2
## Program name does not matter - features gets overwritten
set_io_files(
    {
        chain_id       => $chain_main,
        id             => $id_2,
        file_paths_ref => \@file_paths,
        file_info_href => \%file_info,
    }
);

## Given new program
%io = get_io_files(
    {
        chain_id           => $chain_main,
        id                 => $id,
        file_info_href     => \%file_info,
        order_programs_ref => \@order_programs,
        parameter_href     => \%parameter,
        program_name       => $program_name,
    }
);

## Then infile for new program should be returned
is_deeply( \%{ $file_info{io}{$chain_main}{$id} },
    \%io, q{Got file features for sample_1} );

## Given second program in MAIN chain
my @file_paths_2 = (
    catfile(qw{ a test picard_mergesamfiles file_1.txt}),
    catfile(qw{ a test picard_mergesamfiles file_2.txt}),
);

set_io_files(
    {
        chain_id       => $chain_main,
        id             => $id,
        file_paths_ref => \@file_paths_2,
        file_info_href => \%file_info,
    }
);

## Given the first program in a chain that should inherit
my $first_in_chain_program_name = q{chanjo_sexcheck};
%io = get_io_files(
    {
        chain_id           => $chain_chanjo,
        id                 => $id,
        file_info_href     => \%file_info,
        order_programs_ref => \@order_programs,
        parameter_href     => \%parameter,
        program_name       => $first_in_chain_program_name,
    }
);

is_deeply( \%{ $file_info{io}{$chain_main}{$id} },
    \%io, q{Got inherited file features from MAIN for sample_1} );

## Given a program downstream of PARALLEL chain and other chain
my $downstream_program = q{combinevariantcallsets};

%io = get_io_files(
    {
        chain_id           => $chain_chanjo,
        id                 => $id,
        file_info_href     => \%file_info,
        order_programs_ref => \@order_programs,
        parameter_href     => \%parameter,
        program_name       => $first_in_chain_program_name,
    }
);

is_deeply( \%{ $file_info{io}{$chain_main}{$id} },
    \%io,
    q{Got inherited file features from MAIN from downstream for sample_1} );

## Given a program for sample_2
%io = get_io_files(
    {
        chain_id           => $chain_chanjo,
        id                 => $id_2,
        file_info_href     => \%file_info,
        order_programs_ref => \@order_programs,
        parameter_href     => \%parameter,
        program_name       => $first_in_chain_program_name,
    }
);

## The inherit from sample_2 main chain
is_deeply( \%{ $file_info{io}{$chain_main}{$id_2} },
    \%io, q{Got inherited file features from MAIN for sample_2} );

done_testing();
