#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
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
Readonly my $COMMA                => q{,};
Readonly my $GENOME_BUILD_VERSION => 37;
Readonly my $NEWLINE              => qq{\n};
Readonly my $SPACE                => q{ };

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
    my @modules = (q{MIP::QC::Record});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::QC::Record qw{ add_to_sample_info };

diag(   q{Test add_to_sample_info from Record.pm v}
      . $MIP::QC::Record::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given parameter paths
my %active_parameter = (
    human_genome_reference => catfile(qw{ a test path genome build}),
    log_file               => catfile(qw{ a test dir and log_path}),
    pedigree_file          => catfile(qw{ a test pedigree path }),
);

my %file_info = (
    human_genome_reference_source  => q{GRCh},
    human_genome_reference_version => $GENOME_BUILD_VERSION,
);
my %sample_info;

add_to_sample_info(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        sample_info_href      => \%sample_info,
    }
);

## Then these entries should be set in sample info
is(
    $sample_info{human_genome_build}{path},
    $active_parameter{human_genome_reference},
    q{Added genome build path}
);
is(
    $sample_info{human_genome_build}{source},
    $file_info{human_genome_reference_source},
    q{Added genome build source}
);
is(
    $sample_info{human_genome_build}{version},
    $file_info{human_genome_reference_version},
    q{Added genome build version}
);
is(
    $sample_info{pedigree_file}{path},
    $active_parameter{pedigree_file},
    q{Added pedigree path}
);
is(
    $sample_info{log_file_dir},
    dirname( dirname( $active_parameter{log_file} ) ),
    q{Added log dir path}
);
is(
    $sample_info{last_log_file_path},
    $active_parameter{log_file},
    q{Added log file path}
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
