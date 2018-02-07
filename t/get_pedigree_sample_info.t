#!/usr/bin/env perl

use 5.018;
use Array::Utils qw{ array_diff };
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use List::MoreUtils qw{ any };
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
    my @modules = (q{MIP::Get::Pedigree});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Get::Pedigree qw{ get_sample_info };
use MIP::Log::MIP_log4perl qw{ initiate_logger };

diag(   q{Test get_sample_info from Pedigree.pm v}
      . $MIP::Get::Pedigree::VERSION
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

my %pedigree = (
    family  => q{family_1},
    samples => [
        {
            analysis_type => q{wes},
            father        => 0,
            mother        => 0,
            phenotype     => q{affected},
            sample_id     => q{sample_1},
            sample_origin => q{normal},
            sex           => q{female},
        },
        {
            analysis_type => q{cancer},
            father        => 0,
            mother        => 0,
            phenotype     => q{affected},
            sample_id     => q{sample_2_a},
            sample_origin => q{tumor},
            sex           => q{male},
        },
        {
            analysis_type => q{cancer},
            father        => 0,
            mother        => 0,
            phenotype     => q{unaffected},
            sample_id     => q{sample_2_b},
            sample_origin => q{normal},
            sex           => q{male},
        },
        {
            analysis_type => q{cancer},
            father        => 0,
            mother        => 0,
            phenotype     => q{unaffected},
            sample_id     => q{sample_2_c},
            sample_origin => q{normal},
            sex           => q{male},
        },
        {
            analysis_type => q{wts},
            father        => 0,
            mother        => 0,
            phenotype     => q{unknown},
            sample_id     => q{sample_3},
            sex           => q{other},
        },
        {
            analysis_type => q{wts},
            father        => q{sample_1},
            mother        => q{sample_2},
            phenotype     => q{unknown},
            sample_id     => q{sample_4},
            sex           => q{unknown},
        },
    ],
);

## Test for a single sample entry
my @sample_info = get_sample_info(
    {
        get_values_for_key          => q{sample_id},
        pedigree_href               => \%pedigree,
        sample_info_intersect_key   => q{sample_origin},
        sample_info_intersect_value => q{tumor},
    }
);

my @output_info = (qw{sample_2_a});

is( array_diff( @output_info, @sample_info ),
    0, q{Test for a single info retrieval} );

## Test for multiple sample entry
@sample_info = get_sample_info(
    {
        get_values_for_key          => q{sample_id},
        pedigree_href               => \%pedigree,
        sample_info_intersect_key   => q{analysis_type},
        sample_info_intersect_value => q{cancer},
    }
);

@output_info = (qw{sample_2_a sample_2_b sample_2_c});

is( array_diff( @output_info, @sample_info ),
    0, q{Test for multiple info retrieval} );

## Test for empty output
@sample_info = get_sample_info(
    {
        get_values_for_key          => q{sample_origin},
        pedigree_href               => \%pedigree,
        sample_info_intersect_key   => q{analysis_type},
        sample_info_intersect_value => q{wts},
    }
);

@output_info = ();

is( array_diff( @output_info, @sample_info ),
    0, q{Test for non-existing entry retrieval} );

## Test for default value of get_values_for_key
@sample_info = get_sample_info(
    {
        pedigree_href               => \%pedigree,
        sample_info_intersect_key   => q{analysis_type},
        sample_info_intersect_value => q{cancer},
    }
);

@output_info = (qw{sample_2_a sample_2_b sample_2_c});

is( array_diff( @output_info, @sample_info ),
    0, q{Test for default output for get_values_for_key} );

## Test for default output, i.e. all sample_ids
@sample_info = get_sample_info(
    {
        pedigree_href => \%pedigree,
    }
);

@output_info =
  (qw{ sample_1 sample_2_a sample_2_b sample_2_c sample_3 sample_4 });

is( array_diff( @output_info, @sample_info ),
    0, q{Test for default output for all sample_ids} );

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
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}

