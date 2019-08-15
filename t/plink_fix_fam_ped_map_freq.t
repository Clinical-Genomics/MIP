#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
use Params::Check qw{ check allow last_error };
use Test::More;
use warnings qw{ FATAL utf8 };
use utf8;
use 5.026;

## CPANM
use autodie;
use Modern::Perl qw{ 2018 };
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
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Program::Variantcalling::Plink});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Variantcalling::Plink qw{ plink_fix_fam_ped_map_freq };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test plink_fix_fam_ped_map_freq from Plink.pm v}
      . $MIP::Program::Variantcalling::Plink::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ plink2 };

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    binary_fileset_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--bfile}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    fam_file_path => {
        input           => catfile(qw{ outcase_file_directory case_id .fam }),
        expected_output => q{--fam}
          . $SPACE
          . catfile(qw{ outcase_file_directory case_id .fam }),
    },
    outfile_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--out}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
);

my %specific_argument = (
    allow_no_sex => {
        input           => 1,
        expected_output => q{--allow-no-sex},
    },
    binary_fileset_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--bfile}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    fam_file_path => {
        input           => catfile(qw{ outcase_file_directory case_id .fam }),
        expected_output => q{--fam}
          . $SPACE
          . catfile(qw{ outcase_file_directory case_id .fam }),
    },
    freqx => {
        input           => 1,
        expected_output => q{--freqx},
    },
    make_just_fam => {
        input           => 1,
        expected_output => q{--make-just-fam},
    },
    outfile_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--out}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    recode => {
        input           => 1,
        expected_output => q{--recode},
    },
);

# Coderef - enables generalized use of generate call
my $module_function_cref = \&plink_fix_fam_ped_map_freq;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            required_argument_href     => \%required_argument,
            module_function_cref       => $module_function_cref,
            function_base_commands_ref => \@function_base_commands,
            do_test_base_command       => 1,
        }
    );
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
