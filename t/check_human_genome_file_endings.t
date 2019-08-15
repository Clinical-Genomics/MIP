#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
        q{MIP::Set::Parameter}    => [qw{ set_human_genome_reference_features }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Check::Reference});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Reference qw{ check_human_genome_file_endings };
use MIP::Set::Parameter qw{ set_human_genome_reference_features };

diag(   q{Test check_human_genome_file_endings from Reference.pm v}
      . $MIP::Check::Reference::VERSION
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

my %parameter;

my %active_parameter = (
    human_genome_reference =>
      catfile( $Bin, qw{ data references grch37_homo_sapiens_-d5-.fasta } ),
    reference_dir => catfile( $Bin, qw{ data references } ),
);
## File info hash
my %file_info = (

    # Human genome meta files
    human_genome_reference_file_endings => [qw{ .dict .fai }],
);

## Detect version and source of the human_genome_reference: Source (hg19 or GRCh).
set_human_genome_reference_features(
    {
        file_info_href         => \%file_info,
        human_genome_reference => basename( $active_parameter{human_genome_reference} ),
        log                    => $log,
        parameter_href         => \%parameter,
    }
);

check_human_genome_file_endings(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        log                   => $log,
        parameter_href        => \%parameter,
        parameter_name        => q{human_genome_reference},
    }
);

is( $parameter{human_genome_reference}{build_file},
    0, q{Set build file switch for human genome reference to 0} );

$active_parameter{human_genome_reference} = q{not an existing reference};

check_human_genome_file_endings(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        log                   => $log,
        parameter_href        => \%parameter,
        parameter_name        => q{human_genome_reference},
    }
);

is( $parameter{human_genome_reference}{build_file},
    1, q{Set build file switch for human genome reference to 1} );

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
