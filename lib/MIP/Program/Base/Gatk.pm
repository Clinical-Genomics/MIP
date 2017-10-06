package MIP::Program::Base::Gatk;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gatk_base };
}

## Constants
Readonly my $SPACE => q{ };

sub gatk_base {

## pgatk_base

## Function : Perl wrapper for Gatk base. Based on Gatk v?????
## Returns  : @commands
## Arguments: $commands_ref, $analysis_type, $intervals_ref, $referencefile_path, $pedigree, $pedigree_validation_type, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level, $FILEHANDLE
##          : $commands_ref                          => List of commands added earlier
##          : $analysis_type                         => Analysis type
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $referencefile_path                    => Reference sequence file
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Logging level
##          : $FILEHANDLE                            => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $analysis_type;
    my $intervals_ref;
    my $referencefile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $downsample_to_coverage;
    my $FILEHANDLE;

    ## Default(s)
    my $logging_level;
    my $gatk_disable_auto_index_and_file_lock;

    my $tmpl = {
        commands_ref =>
          { default => [], strict_type => 1, store => \$commands_ref },
        analysis_type => {
            defined     => 1,
            strict_type => 1,
            store       => \$analysis_type
        },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        pedigree                 => { strict_type => 1, store => \$pedigree },
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT}],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        downsample_to_coverage => {
            strict_type => 1,
            store       => \$downsample_to_coverage
        },
        logging_level => {
            default     => q{INFO},
            allow       => [qw{ INFO ERROR FATAL }],
            strict_type => 1,
            store       => \$logging_level
        },
        gatk_disable_auto_index_and_file_lock => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$gatk_disable_auto_index_and_file_lock
        },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = @{$commands_ref};

    if ($analysis_type) {

        push @commands, q{--analysis_type} . $SPACE . $analysis_type;
    }

    push @commands, q{--logging_level} . $SPACE . $logging_level;

# GATK CatVariants does not have an analysis_type parameter and doesn't use pedigree_validation_type
    if ( $pedigree_validation_type && $analysis_type ) {

        push @commands,
          q{--pedigreeValidationType} . $SPACE . $pedigree_validation_type;
    }

    if ($pedigree) {

        push @commands, q{--pedigree} . $SPACE . $pedigree;
    }

    if ($downsample_to_coverage) {

        push @commands,
          q{--downsample_to_coverage} . $SPACE . $downsample_to_coverage;
    }

    if ($gatk_disable_auto_index_and_file_lock) {

        push @commands,
          q{--disable_auto_index_creation_and_locking_when_reading_rods};
    }

    if ( @{$intervals_ref} ) {

        push @commands,
          q{--intervals} . $SPACE . join $SPACE . q{--intervals} . $SPACE,
          @{$intervals_ref};
    }

    if ($referencefile_path) {

      if($analysis_type) {

        push @commands, q{--reference_sequence} . $SPACE . $referencefile_path;
      }
      else {

        push @commands, q{--reference} . $SPACE . $referencefile_path;
      }
    }


    return @commands;
}

1;
