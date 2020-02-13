#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Log::MIP_log4perl qw{ initiate_logger };
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
    my %perl_module = (
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
        q{MIP::Script::Utils}     => [qw{ help }],
    );

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

use MIP::Update::Contigs qw{ size_sort_select_file_contigs };

diag(   q{Test size_sort_select_file_contigs from Contigs.pm v}
      . $MIP::Update::Contigs::VERSION
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

my $consensus_analysis_type = q{wes};

## Given a reference hash of array, when no hash of array to sort exists
my %file_info = ( contigs_size_ordered => [qw{ chr1 chr2 chr3 chrM}], );

trap {
    size_sort_select_file_contigs(
        {
            consensus_analysis_type => $consensus_analysis_type,
            file_info_href          => \%file_info,
            hash_key_sort_reference => q{contigs_size_ordered},
            hash_key_to_sort        => q{select_file_contigs},
            log                     => $log,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if hash key to sort does not exist in supplied hash} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if hash key to sort does not exist in supplied hash} );

## Given a hash of array to sort, when no hash of array as reference exists
%file_info = ( select_file_contigs => [qw{chr2 chr1 chrM chr3}], );

trap {
    size_sort_select_file_contigs(
        {
            consensus_analysis_type => $consensus_analysis_type,
            file_info_href          => \%file_info,
            hash_key_sort_reference => q{contigs_size_ordered},
            hash_key_to_sort        => q{select_file_contigs},
            log                     => $log,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if hash key for reference does not exist in supplied hash} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if hash key for reference does not exist in supplied hash}
);

## Given wes, when lacking mitochondrial contig in array to sort
%file_info = (
    contigs_size_ordered => [qw{ chr1 chr2 chr3 chrM}],
    select_file_contigs  => [qw{chr2 chr1 chr3}],
);

@{ $file_info{sorted_select_file_contigs} } = size_sort_select_file_contigs(
    {
        consensus_analysis_type => $consensus_analysis_type,
        file_info_href          => \%file_info,
        hash_key_sort_reference => q{contigs_size_ordered},
        hash_key_to_sort        => q{select_file_contigs},
        log                     => $log,
    }
);

## Remove chrM for comparison
pop @{ $file_info{contigs_size_ordered} };

## Then return sorted array according to reference (minus chrM in ref)
is_deeply(
    \@{ $file_info{contigs_size_ordered} },
    \@{ $file_info{sorted_select_file_contigs} },
    q{Sorted contigs according to reference for wes}
);

## Given wes, when lacking a contig in reference array
%file_info = (
    contigs_size_ordered => [qw{ chr2 chr3 chrM}],
    select_file_contigs  => [qw{chr2 chr1 chrM chr3}],
);

trap {
    size_sort_select_file_contigs(
        {
            consensus_analysis_type => $consensus_analysis_type,
            file_info_href          => \%file_info,
            hash_key_sort_reference => q{contigs_size_ordered},
            hash_key_to_sort        => q{select_file_contigs},
            log                     => $log,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if contig is lacking from sorted array due to reference} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if contig is lacking from sorted array due to reference} );

## Given all ok parameters
%file_info = (
    contigs_size_ordered => [qw{ chr1 chr2 chr3 chrM}],
    select_file_contigs  => [qw{chr2 chr1 chrM chr3}],
);

@{ $file_info{sorted_select_file_contigs} } = size_sort_select_file_contigs(
    {
        consensus_analysis_type => $consensus_analysis_type,
        file_info_href          => \%file_info,
        hash_key_sort_reference => q{contigs_size_ordered},
        hash_key_to_sort        => q{select_file_contigs},
        log                     => $log,
    }
);

## Then return a complete sorted array
is_deeply(
    \@{ $file_info{contigs_size_ordered} },
    \@{ $file_info{sorted_select_file_contigs} },
    q{Sorted contigs according to reference}
);

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
