#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname  };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = 1.0.2;

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
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Program::Variantcalling::Vep});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Variantcalling::Vep qw{variant_effect_predictor};
use MIP::Test::Commands qw{test_function};

diag(   q{Test variant_effect_predictor from Vep.pm v}
      . $MIP::Program::Variantcalling::Vep::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $VARIANT_BUFFERT_SIZE => 20_000;

## Base arguments
my @function_base_commands = qw{ vep };

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1> stdoutfile.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

my %specific_argument = (
    assembly => {
        input           => q{GRCh37},
        expected_output => q{--assembly} . $SPACE . q{GRCh37},
    },
    buffer_size => {
        input           => $VARIANT_BUFFERT_SIZE,
        expected_output => q{--buffer_size} . $SPACE . $VARIANT_BUFFERT_SIZE,
    },
    cache_directory => {
        input           => catdir( q{test_dir}, q{test_cache_dir} ),
        expected_output => q{--dir_cache}
          . $SPACE
          . catdir( q{test_dir}, q{test_cache_dir} ),
    },
    distance => {
        input           => 10,
        expected_output => q{--distance} . $SPACE . q{10},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    fork => {
        input           => 1,
        expected_output => q{--fork} . $SPACE . q{1},
    },
    infile_format => {
        input           => q{vcf},
        expected_output => q{--format} . $SPACE . q{vcf},
    },
    infile_path => {
        input           => catfile( q{test_dir}, q{infile.vcf} ),
        expected_output => q{--input_file}
          . $SPACE
          . catfile( q{test_dir}, q{infile.vcf} ),
    },
    outfile_format => {
        input           => q{vcf},
        expected_output => q{--} . q{vcf},
    },
    outfile_path => {
        input           => catfile( q{test_dir}, q{infile.vcf} ),
        expected_output => q{--output_file}
          . $SPACE
          . catfile( q{test_dir}, q{infile.vcf} ),
    },
    plugins_dir_path => {
        input           => catdir(qw{ test_dir plugins }),
        expected_output => q{--dir_plugins}
          . $SPACE
          . catdir(qw{ test_dir plugins }),
    },
    plugins_ref => {
        inputs_ref      => [qw{ LoFtool LoF }],
        expected_output => q{--plugin LoFtool} . $SPACE . q{--plugin LoF},
    },
    reference_path => {
        input           => catfile( q{test_dir}, q{hum_ref.pl} ),
        expected_output => q{--fasta}
          . $SPACE
          . catfile( q{test_dir}, q{hum_ref.pl} ),
    },
    regions_ref => {
        inputs_ref      => [qw{ 1 2 }],
        expected_output => q{--chr} . $SPACE . q{1,2},
    },
    vep_features_ref => {
        inputs_ref      => [qw{ tsl hgvs}],
        expected_output => q{--tsl} . $SPACE . q{--hgvs},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&variant_effect_predictor;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

HASHES_OF_ARGUMENTS:
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

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $program_name
##         : $program_name => Name of the script

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
