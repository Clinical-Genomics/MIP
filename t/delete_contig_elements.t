#!/usr/bin/env perl

use 5.026;
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
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

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
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
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
    my @modules = (q{MIP::Delete::List});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Delete::List qw{ delete_contig_elements };

diag(   q{Test delete_contig_elements from List.pm v}
      . $MIP::Delete::List::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given contigs, when no prefix
my @contigs        = qw{ 1 2 3 4 Y};
my @remove_contigs = qw{ 2 4 Y};

my @cleansed_contigs = delete_contig_elements(
    {
        elements_ref       => \@contigs,
        remove_contigs_ref => \@remove_contigs,
    }
);

## Define expected outcome
my @expected_contigs = qw{ 1 3 };

## Then remove the contigs
is_deeply( \@cleansed_contigs, \@expected_contigs, q{Removed contigs} );

## Given contigs, when prefix
my @chr_contigs = qw{ chr1 chr2 chr3 chr4 chrY};

@cleansed_contigs = delete_contig_elements(
    {
        elements_ref       => \@chr_contigs,
        remove_contigs_ref => \@remove_contigs,
    }
);

## Define expected outcome
@expected_contigs = qw{ chr1 chr3 };

## Then remove the contigs irrespective of prefix
is_deeply( \@cleansed_contigs, \@expected_contigs, q{Removed contigs with chr prefix} );

## Given contigs, when prefix in remove
@chr_contigs = qw{ chr1 chr2 chr3 chr4 chrY};
my @chr_remove_contigs = qw{ chr1 chrY };

@cleansed_contigs = delete_contig_elements(
    {
        elements_ref       => \@chr_contigs,
        remove_contigs_ref => \@chr_remove_contigs,
    }
);

## Define expected outcome
@expected_contigs = qw{ chr2 chr3 chr4 };

## Then remove the contigs irrespective of prefix
is_deeply( \@cleansed_contigs, \@expected_contigs,
    q{Removed contigs independent of prefix} );

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
