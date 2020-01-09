#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use Test::More;
use warnings qw{ FATAL utf8 };
use utf8;
use 5.026;

## CPANM
use autodie;
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.2';

## Constants
Readonly my $COMMA               => q{,};
Readonly my $NEWLINE             => qq{\n};
Readonly my $AUTOSOMAL_CONTIG_NR => 22;
Readonly my $SPACE               => q{ };

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
    my %perl_module = (
        q{MIP::Script::Utils}     => [qw{ help }],
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Get::File});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Get::File qw{ get_select_file_contigs };
use MIP::Log::MIP_log4perl qw{ initiate_logger };

diag(   q{Test get_select_file_contigs from File.pm v}
      . $MIP::Get::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $test_dir      = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

## Given proper input data
my %file_info;
my $select_file_path =
  catfile( $Bin, qw{ data 643594-miptest aggregated_gene_panel_test.txt } );

@{ $file_info{select_file_contigs} } = get_select_file_contigs(
    {
        select_file_path => $select_file_path,
        log              => $log,
    }
);
my @expected_contigs = ( 1 .. $AUTOSOMAL_CONTIG_NR, qw{ X Y MT} );

## Then return the expected contigs
is_deeply( \@{ $file_info{select_file_contigs} },
    \@expected_contigs, q{Got select file contigs} );

## Given inproper file path
my $wrong_file = catfile( $Bin, qw{ data 643594-miptest 643594-miptest_pedigree.yaml } );

trap {
    get_select_file_contigs(
        {
            select_file_path => $wrong_file,
            log              => $log,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if contigs cannot be found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if contigs cannot be found} );

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
