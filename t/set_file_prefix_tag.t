#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
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
    my @modules = (q{MIP::Set::File});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Set::File qw{ set_file_prefix_tag };

diag(   q{Test set_file_prefix_tag from File.pm v}
      . $MIP::Set::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $sample_id     = q{sample-1};
my $current_chain = q{MAIN};

my @order_parameters = qw{ bwa_mem pmerge pmark };
my %active_parameter = (
    bwa_mem => 1,
    manta   => 1,
    pmark    => 1,
    pmerge   => 0,
);
my %expected_file_tag;
my %file_info;
my %parameter = (
    bwa_mem => q{mem},
    manta   => q{manta},
    pmark    => q{md},
    pmerge   => q{merge},
);
my %temp_file_ending;

## Given MAIN chain file tags, when merge is turned off
PROGRAM:
foreach my $program (@order_parameters) {

    $temp_file_ending{$current_chain}{$sample_id} = set_file_prefix_tag(
        {
            active_parameter_href => \%active_parameter,
            current_chain         => $current_chain,
            file_tag              => $parameter{$program},
            file_info_href        => \%file_info,
            id                    => $sample_id,
            mip_program_name      => $program,
            temp_file_ending_href => \%temp_file_ending,
        }
    );

}

## Define what to expect
# First file tag
$expected_file_tag{$sample_id}{bwa_mem}{file_tag} = $parameter{bwa_mem};

# Propagated
$expected_file_tag{$sample_id}{pmerge}{file_tag} = $parameter{bwa_mem};

# Sequential
$expected_file_tag{$sample_id}{pmark}{file_tag} =
  $parameter{bwa_mem} . $parameter{pmark};

## Then 3 file tags should be added where one is sequential and one is just propagated
is_deeply( \%file_info, \%expected_file_tag,
    q{Added file prefix tags for MAIN } );

## Given other chain than MAIN
push @order_parameters, q{manta};
my $other_chain = q{SV};

# Clear previous file tag builds
%file_info        = ();
%temp_file_ending = ();

## Given SV chain file tags, when merge is turned off
PROGRAM:
foreach my $program (@order_parameters) {

    $temp_file_ending{$current_chain}{$sample_id} = set_file_prefix_tag(
        {
            active_parameter_href => \%active_parameter,
            current_chain         => $current_chain,
            file_tag              => $parameter{$program},
            file_info_href        => \%file_info,
            id                    => $sample_id,
            mip_program_name      => $program,
            temp_file_ending_href => \%temp_file_ending,
        }
    );
}
$expected_file_tag{$sample_id}{manta}{file_tag} =
  $expected_file_tag{$sample_id}{pmark}{file_tag} . $parameter{manta};

is_deeply( \%file_info, \%expected_file_tag,
    q{Added file prefix tags for chain that inherits from MAIN } );
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
