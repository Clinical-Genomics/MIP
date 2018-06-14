#!/usr/bin/env perl
use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
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
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Program::Variantcalling::BlobFish});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Variantcalling::BlobFish qw{ blobfish_all_vs_all };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test blobfish_all_vs_all from BlobFish.pm v}
      . $MIP::Program::Variantcalling::BlobFish::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $function_base_command = q{BlobFish.py --allvsall};

my %base_argument = (
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1> stdoutfile.test},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => $function_base_command,
    },
);

my %required_argument = (
    outdir_path => {
        input           => catfile(qw{ output directory }),
        expected_output => q{--dir} . $SPACE . catfile(qw{ output directory }),
    },
    tx_to_gene_path => {
        input           => catfile(qw{ path to gene2tx.csv }),
        expected_output => q{--tx}
          . $SPACE
          . catfile(qw{ path to gene2tx.csv }),
    },
    infile_paths_ref => {
        inputs_ref => [
            catfile(qw{ dir A1 quant.sf }), catfile(qw{ dir A2 quant.sf }),
            catfile(qw{ dir B1 quant.sf }), catfile(qw{ dir B2 quant.sf })
        ],
        expected_output => q{--paths}
          . $SPACE
          . catfile(qw{ dir A1 quant.sf })
          . $SPACE
          . catfile(qw{ dir A2 quant.sf })
          . $SPACE
          . catfile(qw{ dir B1 quant.sf })
          . $SPACE
          . catfile(qw{ dir B2 quant.sf }),
    },
    insamples_ref => {
        inputs_ref      => [ q{A}, q{A}, q{B}, q{B} ],
        expected_output => q{--sample}
          . $SPACE . q{B}
          . $SPACE . q{A}
          . $SPACE . q{B}
          . $SPACE . q{B},
    },
);

my %specific_argument = (
    outdir_path => {
        input           => catfile(qw{ output directory }),
        expected_output => q{--dir} . $SPACE . catfile(qw{ output directory }),
    },
    tx_to_gene_path => {
        input           => catfile(qw{ path to gene2tx.csv }),
        expected_output => q{--tx}
          . $SPACE
          . catfile(qw{ path to gene2tx.csv }),
    },
    infile_paths_ref => {
        inputs_ref => [
            catfile(qw{ dir A1 quant.sf }), catfile(qw{ dir A2 quant.sf }),
            catfile(qw{ dir B1 quant.sf }), catfile(qw{ dir B2 quant.sf })
        ],
        expected_output => q{--paths}
          . $SPACE
          . catfile(qw{ dir A1 quant.sf })
          . $SPACE
          . catfile(qw{ dir A2 quant.sf })
          . $SPACE
          . catfile(qw{ dir B1 quant.sf })
          . $SPACE
          . catfile(qw{ dir B2 quant.sf }),
    },
    insamples_ref => {
        inputs_ref      => [ q{A}, q{A}, q{B}, q{B} ],
        expected_output => q{--sample}
          . $SPACE . q{A}
          . $SPACE . q{A}
          . $SPACE . q{B}
          . $SPACE . q{B},
    },
);

my $module_function_cref = \&blobfish_all_vs_all;

my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href          => $argument_href,
            required_argument_href => \%required_argument,
            module_function_cref   => $module_function_cref,
            function_base_command  => $function_base_command,
            do_test_base_command   => 1,
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
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
