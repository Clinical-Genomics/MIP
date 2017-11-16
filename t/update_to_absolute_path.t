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
use utf8;
use warnings qw{ FATAL utf8 };
use 5.018;

## CPANM
use autodie;
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

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
    my @modules = (q{MIP::Update::Path});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Update::Path qw{ update_to_absolute_path };

diag(   q{Test update_to_absolute_path from Update::Path.pm v}
      . $MIP::Update::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %parameter =
  ( dynamic_parameter =>
      { absolute_path => [qw{ hash array scalar no_absolute_path }] } );

my %active_parameter = (
    hash => {
        scalar =>
          catfile(qw{ data reference GRCh37_homo_sapiens_-d5-.fasta.gz })
    },
    array =>
      [ catfile(qw{ data reference GRCh37_homo_sapiens_-d5-.fasta.gz }) ],
    scalar => catfile(qw{ data reference GRCh37_homo_sapiens_-d5-.fasta.gz }),
    no_absolute_path => catfile(qw{ not a file }),
);

update_to_absolute_path(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
    }
);

## NOTE: Update_to_absolute_path uppdated path for hash key and not value
my $hash_key_path = catfile( $Bin, qw{ scalar } );
my $expected_value =
  catfile(qw{ data reference GRCh37_homo_sapiens_-d5-.fasta.gz });
is( $active_parameter{hash}{$hash_key_path},
    $expected_value, q{Set hash absolute path} );

is( $active_parameter{array}[0], $expected_value, q{Set array absolute path} );

is( $active_parameter{scalar}, $expected_value, q{Set scalar absolute path} );

my %no_path = ( no_absolute_path => catfile(qw{ not a file }), );

update_to_absolute_path(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%no_path,
    }
);

isnt( $no_path{no_absolute_path}, undef, q{Got error message} );

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
