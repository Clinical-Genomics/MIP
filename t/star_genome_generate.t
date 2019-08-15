#!/usr/bin/env perl

use 5.026;
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
    my @modules = (q{MIP::Program::Alignment::Star});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Alignment::Star qw{ star_genome_generate };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test star_genome_generate from Star.pm v}
      . $MIP::Program::Alignment::Star::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $READ_LENGTH   => 150;
Readonly my $THREAD_NUMBER => 16;

my @function_base_commands = qw{ STAR --runMode genomeGenerate };

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
        expected_output => \@function_base_commands,
    },
);

my %required_argument = (
    fasta_path => {
        input           => catfile(qw{ dir test_file.fasta }),
        expected_output => q{--genomeFastaFiles}
          . $SPACE
          . catfile(qw{ dir test_file.fasta }),
    },
    genome_dir_path => {
        input           => catfile(qw{ dir genome_dir_path }),
        expected_output => q{--genomeDir} . $SPACE . catfile(qw{ dir genome_dir_path }),
    },
    gtf_path => {
        input           => catfile(qw{ dir test_gtf.gtf }),
        expected_output => q{--sjdbGTFfile} . $SPACE . catfile(qw{ dir test_gtf.gtf }),
    },
);

my %specific_argument = (
    fasta_path => {
        input           => catfile(qw{ dir test_file.fasta }),
        expected_output => q{--genomeFastaFiles}
          . $SPACE
          . catfile(qw{ dir test_file.fasta }),
    },
    genome_dir_path => {
        input           => catfile(qw{ dir genome_dir_path }),
        expected_output => q{--genomeDir} . $SPACE . catfile(qw{ dir genome_dir_path }),
    },
    gtf_path => {
        input           => catfile(qw{ dir test_gtf.gtf }),
        expected_output => q{--sjdbGTFfile} . $SPACE . catfile(qw{ dir test_gtf.gtf }),
    },
    read_length => {
        input           => $READ_LENGTH,
        expected_output => q{--sjdbOverhang} . $SPACE . $READ_LENGTH,
    },
    thread_number => {
        input           => $THREAD_NUMBER,
        expected_output => q{--runThreadN} . $SPACE . $THREAD_NUMBER,
    },
);

my $module_function_cref = \&star_genome_generate;

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
