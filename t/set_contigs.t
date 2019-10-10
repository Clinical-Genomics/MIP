#!/usr/bin/env perl

use Modern::Perl qw{ 2014 };
use warnings qw{ FATAL utf8 };
use autodie;
use 5.026;
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir };
use Getopt::Long;
use Test::More;

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.1';

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};
Readonly my $COMMA   => q{,};

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
    my %perl_module;

    $perl_module{q{MIP::Script::Utils}} = [qw{ help }];

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Set::Contigs});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Set::Contigs qw{ set_contigs };

diag(   q{Test set_contigs from List.pm v}
      . $MIP::Set::Contigs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $INDEX_SIZE_ORDERED_CHR_X => 7;

my %file_info;

my @refseq_contigs = qw{
  chr1 chr2 chr3 chr4 chr5 chr6
  chr7 chr8 chr9 chr10 chr11 chr12
  chr13 chr14 chr15 chr16 chr17 chr18
  chr19 chr20 chr21 chr22 chrX chrY
  chrM };

my @ensembl_contigs = qw{
  1 2 3 4 5 6 7 8 9 10
  11 12 13 14 15 16 17 18 19 20
  21 22 X Y MT };

my @grch37_contigs = qw{
  1 2 3 4 5 6 7 8 9 10
  11 12 13 14 15 16 17 18 19 20
  21 22 X Y MT };
my @grch37_contigs_size_ordered = qw{
  1 2 3 4 5 6 7 X 8 9
  10 11 12 13 14 15 16 17 18 19
  20 21 22 Y MT };

## Tests

# grch
set_contigs(
    {
        file_info_href         => \%file_info,
        human_genome_reference => q{grch37_homo_sapiens_-d5-.fasta},
    }
);

is( $file_info{contigs}[-1], q{MT}, q{Set grch reference contigs} );

is( $file_info{contigs_size_ordered}[$INDEX_SIZE_ORDERED_CHR_X],
    q{X}, q{Set grch reference size ordered contigs} );

is_deeply( \@{ $file_info{bam_contigs} },
    \@grch37_contigs, q{Set grch37 reference bam contigs} );
is_deeply(
    \@{ $file_info{bam_contigs_size_ordered} },
    \@grch37_contigs_size_ordered,
    q{Set grch37 reference size ordered bam contigs}
);

# Hg38
set_contigs(
    {
        file_info_href         => \%file_info,
        human_genome_reference => q{hg38_homo_sapiens_-decoy_hla-.fasta},
    }
);

is( $file_info{contigs}[-1], q{chrM}, q{Set hg38 reference contigs} );

is( $file_info{contigs_size_ordered}[$INDEX_SIZE_ORDERED_CHR_X],
    q{chrX}, q{Set hg38 reference size ordered contigs} );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## build_usage

## Function  : Build the USAGE instructions
## Returns   : ""
## Arguments : $program_name
##          : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
