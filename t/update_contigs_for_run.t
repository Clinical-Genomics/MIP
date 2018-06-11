#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
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
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA         => q{,};
Readonly my $MIN_CONTIG_NR => 1;
Readonly my $MAX_CONTIG_NR => 6;
Readonly my $NEWLINE       => qq{\n};
Readonly my $SPACE         => q{ };

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Update::Contigs});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
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

## Geiven some contigs to remove
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

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help     Display this help message
    -v/--version  Display version
END_USAGE
}
