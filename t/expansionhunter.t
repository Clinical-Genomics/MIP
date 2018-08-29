#!/usr/bin/env perl

use 5.022;
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
our $VERSION = 1.0.0;

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
    my @modules = (q{MIP::Program::Variantcalling::Expansionhunter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Variantcalling::Expansionhunter qw{ expansionhunter };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test expansionhunter from Expansionhunter.pm v}
      . $MIP::Program::Variantcalling::Expansionhunter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Module test constants
Readonly my $MIN_ANCHOR_MAPQ         => 60;
Readonly my $MIN_BASEQ               => 30;
Readonly my $MIN_SCORE               => 0.8;
Readonly my $READ_DEPTH              => 30;
Readonly my $REGION_EXTENSION_LENGTH => 1500;

## Base arguments
my @function_base_commands = qw{ ExpansionHunter };

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
    infile_path => {
        input           => catfile(qw{ a test path }),
        expected_output => q{--bam} . $SPACE . catfile(qw{ a test path }),
    },
    json_outfile_path => {
        input           => catfile(qw{ a test json }),
        expected_output => q{--json} . $SPACE . catfile(qw{ a test json }),
    },
    log_outfile_path => {
        input           => catfile(qw{ a test log }),
        expected_output => q{--log} . $SPACE . catfile(qw{ a test log }),
    },
    reference_genome_path => {
        input           => catfile(qw{ a test fasta }),
        expected_output => q{--ref-fasta}
          . $SPACE
          . catfile(qw{ a test fasta }),
    },
    repeat_specs_dir_path => {
        input           => catdir(qw{ a test rep_spec_dir }),
        expected_output => q{--repeat-specs}
          . $SPACE
          . catdir(qw{ a test rep_spec_dir }),
    },
    vcf_outfile_path => {
        input           => catfile(qw{ a test vcf }),
        expected_output => q{--vcf} . $SPACE . catfile(qw{ a test vcf }),
    },
);

my %specific_argument = (
    infile_path => {
        input           => catfile(qw{ a test path }),
        expected_output => q{--bam} . $SPACE . catfile(qw{ a test path }),
    },
    json_outfile_path => {
        input           => catfile(qw{ a test json }),
        expected_output => q{--json} . $SPACE . catfile(qw{ a test json }),
    },
    log_outfile_path => {
        input           => catfile(qw{ a test log }),
        expected_output => q{--log} . $SPACE . catfile(qw{ a test log }),
    },
    min_anchor_mapq => {
        input           => $MIN_ANCHOR_MAPQ,
        expected_output => q{--min-anchor-mapq} . $SPACE . $MIN_ANCHOR_MAPQ,
    },
    min_baseq => {
        input           => $MIN_BASEQ,
        expected_output => q{--min-baseq} . $SPACE . $MIN_BASEQ,
    },
    min_score => {
        input           => $MIN_SCORE,
        expected_output => q{--min-score} . $SPACE . $MIN_SCORE,
    },
    read_depth => {
        input           => $READ_DEPTH,
        expected_output => q{--read-depth} . $SPACE . $READ_DEPTH,
    },
    reference_genome_path => {
        input           => catfile(qw{ a test fasta }),
        expected_output => q{--ref-fasta}
          . $SPACE
          . catfile(qw{ a test fasta }),
    },
    region_extension_length => {
        input           => $REGION_EXTENSION_LENGTH,
        expected_output => q{--region-extension-length}
          . $SPACE
          . $REGION_EXTENSION_LENGTH,
    },
    repeat_specs_dir_path => {
        input           => catdir(qw{ a test rep_spec_dir }),
        expected_output => q{--repeat-specs}
          . $SPACE
          . catdir(qw{ a test rep_spec_dir }),
    },
    sex => {
        input           => q{female},
        expected_output => q{--sex female},
    },
    vcf_outfile_path => {
        input           => catfile(qw{ a test vcf }),
        expected_output => q{--vcf} . $SPACE . catfile(qw{ a test vcf }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&expansionhunter;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
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
