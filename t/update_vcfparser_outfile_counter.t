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
    my @modules = (q{MIP::Update::Parameters});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Update::Parameters qw{ update_vcfparser_outfile_counter };

diag(   q{Test update_vcfparser_outfile_counter from Parameters.pm v}
      . $MIP::Update::Parameters::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Test the number of oufiles when vcfparser is used with select file and sv_vcfparser without.
my %active_parameter_test = (
    sv_vcfparser          => { type => q{program} },
    vcfparser             => { type => q{program} },
    vcfparser_select_file => 1,
);

update_vcfparser_outfile_counter(
    { active_parameter_href => \%active_parameter_test, } );

is( $active_parameter_test{vcfparser_outfile_count},
    2, q{vcfparser used with a select file -> 2 outfiles. Test passed.} );
is( $active_parameter_test{sv_vcfparser_outfile_count},
    1, q{sv_vcfparser used without a select file -> 1 outfile. Test passed.} );

# Test the number of oufiles when both vcfparser and sv_vcfparser are used with select files.
%active_parameter_test = (
    sv_vcfparser             => { type => q{program} },
    vcfparser                => { type => q{program} },
    sv_vcfparser_select_file => 1,
    vcfparser_select_file    => 1,
);

update_vcfparser_outfile_counter(
    { active_parameter_href => \%active_parameter_test, } );

is( $active_parameter_test{vcfparser_outfile_count},
    2, q{vcfparser used with a select file -> 2 outfiles. Test passed.} );
is( $active_parameter_test{sv_vcfparser_outfile_count},
    2, q{sv_vcfparser used with a select file -> 2 outfiles. Test passed.} );

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
