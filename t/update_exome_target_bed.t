#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
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
    my @modules = (q{MIP::Update::Parameters});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Update::Parameters qw{ update_exome_target_bed };

diag(   q{Test update_exome_target_bed from Parameters.pm v}
      . $MIP::Update::Parameters::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Test hashes
my %file_info = (
    human_genome_reference_version => q{37},
    human_genome_reference_source  => q{grch},
);

my %active_parameter = (
    exome_target_bed => {
        q{genome_reference_source_version_agilent_sureselect_targets_cre_-v1-.bed} =>
          q{sample_1},
        q{test_capture.bed} => q{sample_2},
    },
);

update_exome_target_bed(
    {
        exome_target_bed_file_href     => $active_parameter{exome_target_bed},
        human_genome_reference_source  => $file_info{human_genome_reference_source},
        human_genome_reference_version => $file_info{human_genome_reference_version},
    }
);

foreach my $capture_file ( keys %{ $active_parameter{exome_target_bed} } ) {

    if ( $capture_file ne q{test_capture.bed} ) {

        like( $capture_file, qr/grch/xsm, q{Updated genome source} );

        like( $capture_file, qr/37/xsm, q{Updated genome version} );
    }
    else {

        unlike( $capture_file, qr/grch/xsm, q{Did not updated genome source} );

        unlike( $capture_file, qr/37/xsm, q{Did not updated genome version} );
    }
}

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
