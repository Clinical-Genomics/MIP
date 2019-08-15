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
Readonly my $COMMA         => q{,};
Readonly my $MIN_CONTIG_NR => 1;
Readonly my $MAX_CONTIG_NR => 6;
Readonly my $NEWLINE       => qq{\n};
Readonly my $SPACE         => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Update::Contigs} => [qw{ update_contigs_for_run }],
        q{MIP::Test::Fixtures}  => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Update::Contigs qw{ update_contigs_for_run };

diag(   q{Test update_contigs_for_run from Contigs.pm v}
      . $MIP::Update::Contigs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given some contigs to remove
my $found_male      = 0;
my @exclude_contigs = qw{ 1 2 3 4};

## Define analysis type
my %active_parameter = ( analysis_type => { sample_1 => q{wes}, }, );

my %file_info = (
    contigs_size_ordered => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
    contigs              => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
    select_file_contigs  => [ ( $MIN_CONTIG_NR .. $MAX_CONTIG_NR ), qw{ MT Y} ],
);

update_contigs_for_run(
    {
        file_info_href      => \%file_info,
        exclude_contigs_ref => \@exclude_contigs,
        analysis_type_href  => \%{ $active_parameter{analysis_type} },
        found_male          => $found_male,
        log                 => $log,
    }
);

## Define expected outcome
my %expected_result = (
    contigs_size_ordered => [qw{ 5 6}],
    contigs              => [qw{ 5 6}],
    select_file_contigs  => [qw{ 5 6 }],
);
## Then only keep contig 5 and 6
is_deeply( \%file_info, \%expected_result, q{Updated contigs} );

done_testing();
