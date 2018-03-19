#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use Params::Check qw{ check allow last_error };
use Test::More;
use warnings qw{ FATAL utf8 };
use utf8;
use 5.018;

## CPANM
use autodie;
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $SPACE        => q{ };
Readonly my $NEWLINE      => qq{\n};
Readonly my $COMMA        => q{,};
Readonly my $HG_VERSION   => 19;
Readonly my $GRCH_VERSION => 37;

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
    my @modules = (q{MIP::Set::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Set::Parameter qw{ set_human_genome_reference_features };
use MIP::Log::MIP_log4perl qw{ initiate_logger };

diag(   q{Test set_human_genome_reference_features from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

## Test hash
my %file_info;

## Test Ensemble genome
my $human_genome_reference = q{GRCh37_homo_sapiens_-d5-.fasta};
set_human_genome_reference_features(
    {
        file_info_href         => \%file_info,
        human_genome_reference => $human_genome_reference,
        log                    => $log,
    }
);
is( $file_info{human_genome_reference_version},
    $GRCH_VERSION, q{GRCh version test} );
is( $file_info{human_genome_reference_source}, q{GRCh}, q{GRCh source test} );
is( $file_info{human_genome_reference_name_prefix},
    q{GRCh37_homo_sapiens_-d5-}, q{GRCh prefix test} );
is( $file_info{human_genome_compressed}, 0, q{GRCh compressed test} );

## Test Refseq genome
$human_genome_reference = q{hg19_homo_sapiens.fasta.gz};
set_human_genome_reference_features(
    {
        file_info_href         => \%file_info,
        human_genome_reference => $human_genome_reference,
        log                    => $log,
    }
);
is( $file_info{human_genome_reference_version},
    $HG_VERSION, q{hg version test} );
is( $file_info{human_genome_reference_source}, q{hg}, q{hg source test} );
is( $file_info{human_genome_reference_name_prefix},
    q{hg19_homo_sapiens}, q{hg prefix test} );
is( $file_info{human_genome_compressed}, 1, q{hg compressed test} );

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
