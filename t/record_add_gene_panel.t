#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname  };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
use Log::Log4perl;
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
    my @modules = (q{MIP::QC::Record});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::QC::Record qw{ add_gene_panel };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test add_gene_panel from Record.pm v}
      . $MIP::QC::Record::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $aggregate_gene_panel_file =
  catfile( $Bin, qw{ data 643594-miptest aggregated_gene_panel_test.txt } );
my $aggregate_gene_panels_key = q{select_file};
my $gene_panel                = q{TEST};
my $family_id_test            = q{family_id};
my $program_name_test         = q{vcfparser};
my %sample_info;

my %header_info = (
    display_name => q{gene_panel_test},
    gene_panel   => $gene_panel,
    updated_at   => q{2016-12-08},
    version      => q{1.0},
);

add_gene_panel(
    {
        aggregate_gene_panel_file => $aggregate_gene_panel_file,
        aggregate_gene_panels_key => $aggregate_gene_panels_key,
        family_id                 => $family_id_test,
        program_name              => $program_name_test,
        sample_info_href          => \%sample_info,
    }
);

is(
    exists $sample_info{$program_name_test}{$aggregate_gene_panels_key}
      {gene_panel}{$gene_panel},
    1,
    q{Gene panel key added to $sample_info}
);

while ( my ( $key, $value ) = each %header_info ) {

## Test gene panel info
    my $set_header_value =
      $sample_info{$program_name_test}{$aggregate_gene_panels_key}{gene_panel}
      {$gene_panel}{$key};

    is( $set_header_value, $value,
            q{Gene panel header info value for key: }
          . $key
          . q{ added to $sample_info} );
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
