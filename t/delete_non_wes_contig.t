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
use Modern::Perl qw{ 2014 };
use Readonly;

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
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Delete::List}   => [qw{ delete_non_wes_contig }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Delete::List qw{ delete_non_wes_contig };

diag(   q{Test delete_non_wes_contig from List.pm v}
      . $MIP::Delete::List::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given contigs, when no prefix
my @contigs = qw{ 1 2 3 4 MT};

## Define analysis type
my %active_parameter = ( analysis_type => { sample_1 => q{wes}, }, );

my @no_wes_contigs = delete_non_wes_contig(
    {
        analysis_type_href => \%{ $active_parameter{analysis_type} },
        contigs_ref        => \@contigs,
        contig_names_ref   => [qw{ M MT }],
        log                => $log,
    }
);

## Define expected outcome
my @expected_no_wes_contigs = qw{ 1 2 3 4 };

## Then remove the non wes contigs
is_deeply( \@no_wes_contigs, \@expected_no_wes_contigs, q{Removed non wes contigs} );

## Given contigs, when a non consensus type of run
$active_parameter{analysis_type}{sample_2} = q{wgs};

@no_wes_contigs = delete_non_wes_contig(
    {
        analysis_type_href => \%{ $active_parameter{analysis_type} },
        contigs_ref        => \@contigs,
        log                => $log,
    }
);

## Then remove the non wes contigs
is_deeply( \@no_wes_contigs, \@expected_no_wes_contigs,
    q{Remove non wes contigs for mixed analysis run} );

## Given contigs, when a non consensus type of run
$active_parameter{analysis_type}{sample_1} = q{wgs};

my @has_wes_contigs = delete_non_wes_contig(
    {
        analysis_type_href => \%{ $active_parameter{analysis_type} },
        contigs_ref        => \@contigs,
        log                => $log,
    }
);

## Define expected outcome
my @expected_wes_contigs = qw{ 1 2 3 4 MT };

## Then remove the non wes contigs
is_deeply( \@has_wes_contigs, \@expected_wes_contigs, q{Did not remove non wes contigs} );

done_testing();
