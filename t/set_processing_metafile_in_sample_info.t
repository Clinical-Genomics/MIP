#!/usr/bin/env perl

use Modern::Perl qw{ 2018 };
use warnings qw{FATAL utf8};
use autodie;
use 5.026;    # Require at least perl 5.18
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{-no_match_vars};
use Params::Check qw{check allow last_error};

use FindBin qw{$Bin};    # Find directory of script
use File::Basename qw{dirname basename};
use File::Spec::Functions qw{catfile catdir devnull};
use Getopt::Long;
use Test::More;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{help};

our $USAGE = build_usage( {} );

my $VERBOSE = 0;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

###User Options
GetOptions(
    'h|help' => sub {
        done_testing();
        print {*STDOUT} $USAGE, "\n";
        exit;
    },    #Display help text
    'v|version' => sub {
        done_testing();
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{  } . $VERSION, "\n\n";
        exit;
    },    #Display version number
    'vb|verbose' => $VERBOSE,
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
    my %perl_module = ( q{MIP::Script::Utils} => [qw{help}], );

  PERL_MODULES:
    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

    ## Modules
    my @modules = (qw{MIP::Sample_info});

  MODULES:
    for my $module (@modules) {

        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Sample_info qw{set_processing_metafile_in_sample_info};

diag(   q{Test set_processing_metafile_in_sample_info from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Init hash
my %sample_info;

# Test variables
my $metafile_tag = q{most_complete_bam};
my $metafile     = q{test_metafile};
my $directory    = q{test_directory};
my $path         = catfile( $directory, $metafile );

## Family level
set_processing_metafile_in_sample_info(
    {
        sample_info_href => \%sample_info,
        metafile_tag     => $metafile_tag,
        path             => $path,
    }
);

## Test
is( exists $sample_info{$metafile_tag}, 1, q{Created case level hash key} );

is( $sample_info{$metafile_tag}{path},
    $path, q{Assigned correct value to case level path} );

## Sample level
my $sample_id = q{test_sample_id};

set_processing_metafile_in_sample_info(
    {
        sample_info_href => \%sample_info,
        sample_id        => $sample_id,
        metafile_tag     => $metafile_tag,
        path             => $path,
    }
);

## Test
is( exists $sample_info{sample}{$sample_id}{$metafile_tag},
    1, q{Created sample level hash key} );

is( $sample_info{sample}{$sample_id}{$metafile_tag}{path},
    $path, q{Assigned correct value to sample level path} );

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
