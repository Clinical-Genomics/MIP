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
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::File_info}      => [qw{ set_bam_contigs }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_bam_contigs };

diag(   q{Test set_bam_contigs from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given bam contigs
my $contig_set     = q{bam_contigs};
my @grch38_contigs = qw{
  chr1 chr2 chr3 chr4 chr5 chr6
  chr7 chr8 chr9 chr10 chr11 chr12
  chr13 chr14 chr15 chr16 chr17 chr18
  chr19 chr20 chr21 chr22 chrX chrY
  chrM };

my %file_info;

set_bam_contigs(
    {
        file_info_href      => \%file_info,
        primary_contigs_ref => \@grch38_contigs,
        bam_contig_set_name => $contig_set,
    }
);

## Then
is_deeply( \@{ $file_info{bam_contigs} },
    \@grch38_contigs, q{Set grch38 reference bam contigs} );

done_testing();
