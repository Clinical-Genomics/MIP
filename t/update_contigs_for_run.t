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
use MIP::Test::Fixtures qw{ test_log };

## Constants
Readonly my $MIN_CONTIG_NR => 1;
Readonly my $MAX_CONTIG_NR => 6;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Contigs}        => [qw{ update_contigs_for_run }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Contigs qw{ update_contigs_for_run };

diag(   q{Test update_contigs_for_run from Contigs.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given some contigs to remove
my $include_y       = 0;
my @exclude_contigs = qw{ 1 2 3 4};

my %file_info = (
    bam_contigs              => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
    bam_contigs_size_ordered => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
    contigs                  => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
    contigs_size_ordered     => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
    primary_contigs          => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
    select_file_contigs      => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
);

update_contigs_for_run(
    {
        consensus_analysis_type => q{wes},
        exclude_contigs_ref     => \@exclude_contigs,
        file_info_href          => \%file_info,
        include_y               => $include_y,
    }
);

## Define expected outcome
my %expected_result = (
    bam_contigs              => [qw{ 5 6 Y }],
    bam_contigs_size_ordered => [qw{ 5 6 Y }],
    contigs                  => [qw{ 5 6 }],
    contigs_size_ordered     => [qw{ 5 6 }],
    primary_contigs          => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
    select_file_contigs      => [qw{ 5 6 }],
);

## Then only keep contig 5 and 6 for all except bam contigs
is_deeply( \%file_info, \%expected_result, q{Updated contigs} );

done_testing();
