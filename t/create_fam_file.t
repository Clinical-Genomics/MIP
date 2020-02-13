#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE $TAB };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

my $VERBOSE = 1;
our $VERSION = 1.02;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Pedigree}       => [qw{ create_fam_file }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ create_fam_file };

diag(   q{Test create_fam_file from Pedigree.pm v}
      . $MIP::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create test hashes
my %sample_info_test_hash = (
    sample => {
        q{118-1-2A} => {
            mother    => q{118-2-2U},
            father    => q{118-2-1U},
            sex       => q{female},
            phenotype => q{affected},
            sample_id => q{118-1-2A},
        },
        q{118-2-2U} => {
            mother    => 0,
            father    => 0,
            sex       => q{female},
            phenotype => q{unaffected},
            sample_id => q{118-2-2U},
        },
        q{118-2-1U} => {
            mother    => 0,
            father    => 0,
            sex       => q{male},
            phenotype => q{unaffected},
            sample_id => q{118-2-1U},
        },

    },
    case => q{118},
);

my %active_parameter_test_hash = (
    case_id    => q{118},
    sample_ids => [qw{ 118-1-2A 118-2-2U 118-2-1U }],
  ),

## Create temporary directories and paths
  my $test_dir = File::Temp->newdir();
my $test_fam_file_path = catfile( $test_dir, q{test_file.fam} );
my $test_sbatch_path   = catfile( $test_dir, q{test_file.sbatch} );
my $filehandle;

## Create temp logger for Pedigree.pm
my $test_log_path = catfile( $test_dir, q{test.log} );
$active_parameter_test_hash{log_file} = $test_log_path;

my $log = test_log( {} );

my @execution_modes = qw{ system sbatch };
for my $execution_mode (@execution_modes) {

    if ( $execution_mode eq q{sbatch} ) {

        $filehandle = IO::Handle->new();
        open $filehandle, q{>}, $test_sbatch_path
          or croak q{Could not open filehandle};

        # Build shebang
        my @commands = (q{#!/usr/bin/env bash});
        unix_write_to_file(
            {
                commands_ref => \@commands,
                filehandle   => $filehandle,
                separator    => $SPACE,
            }
        );
    }

    # Run the create fam file test
    create_fam_file(
        {
            active_parameter_href => \%active_parameter_test_hash,
            execution_mode        => $execution_mode,
            fam_file_path         => $test_fam_file_path,
            filehandle            => $filehandle,
            log                   => $log,
            sample_info_href      => \%sample_info_test_hash,
        }
    );

    if ( $execution_mode eq q{sbatch} ) {

        ## Then fam sbatch file should have been created
        ok( -e $test_sbatch_path, q{fam command written to sbatch file} );

        # Execute the sbatch
        close $filehandle;
        system qq{bash $test_sbatch_path};
        unlink $test_sbatch_path or carp qq{Could not unlink $test_sbatch_path};
    }

    ## Then fam file should have been created
    ok( -e $test_fam_file_path, q{fam file created} );

    # Given fam headers
    open my $fh, q{<}, $test_fam_file_path
      or croak qq{Could not open $test_fam_file_path for reading};

    ## Get headers
    my $file_header = <$fh>;
    chomp $file_header;
    my $expected_header = qq{#family_id\tsample_id\tfather\tmother\tsex\tphenotype};

    ## Then header should be found in file
    is( $file_header, $expected_header, q{header included in fam} );

    ## Testing pattern of created pedigree file

    # Creating array of expected pedigree lines
    my @expected_pedigree_lines;

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_test_hash{sample_ids} } ) {

        my $sample_line = $active_parameter_test_hash{case_id};

      HEADER:
        foreach my $header ( split $TAB, $expected_header ) {

            if ( defined $sample_info_test_hash{sample}{$sample_id}{$header} ) {
                $sample_line .=
                  $TAB . $sample_info_test_hash{sample}{$sample_id}{$header};
            }
        }
        push @expected_pedigree_lines, $sample_line;
    }

    # Reading pedigree lines from fam file into array for testing
    my @pedigree_lines = <$fh>;
    close $fh;
    chomp @pedigree_lines;

    ## Then correct pedigree lines should be found
    is( @pedigree_lines, @expected_pedigree_lines, q{fam file has correct information} );

    # If the fam file exists when running in sbatch mode no sbatch file will be created
    unlink $test_fam_file_path
      or carp qq{Could not unlink $test_fam_file_path};
}

## Given no filehandle when in sbatch mode
# Run the create fam file test
trap {
    create_fam_file(
        {
            active_parameter_href => \%active_parameter_test_hash,
            execution_mode        => q{sbatch},
            fam_file_path         => $test_fam_file_path,
            log                   => $log,
            sample_info_href      => \%sample_info_test_hash,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if no filehandle supplied in sbatch mode} );
like(
    $trap->stderr,
    qr/Please \s+ supply \s+ filehandle \s+ to/xms,
    q{Throw fatal log message if no filehandle supplied in sbatch mode}
);

done_testing();
