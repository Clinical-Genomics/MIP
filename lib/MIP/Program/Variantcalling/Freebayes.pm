package MIP::Program::Variantcalling::Freebayes;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ freebayes_calling };
}

## Constants
Readonly my $SPACE => q{ };

sub freebayes_calling {

## Function : Perl wrapper for generic commands module.
## Returns  : @commands

## Arguments: $apply_standard_filter      => Use stringent input base and mapping quality filters. Equivalent to -m 30 -q 20 -R 0 -S 0
##          : $calculate_genotype_quality => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output
##          : $FILEHANDLE                 => Filehandle to write to
##          : $infile_paths_ref           => Infile paths {REF}
##          : $referencefile_path         => Reference sequence file
##          : $stdoutfile_path            => Stdoutfile path
##          : $stderrfile_path            => Stderrfile path
##          : $use_best_n_alleles         => Reduce the number of alleles that are considered,
    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $apply_standard_filter;
    my $calculate_genotype_quality;
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $referencefile_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $use_best_n_alleles;

    my $tmpl = {
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        FILEHANDLE      => { store => \$FILEHANDLE, },
        apply_standard_filter => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$apply_standard_filter,
            strict_type => 1,
        },
        calculate_genotype_quality => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$calculate_genotype_quality,
            strict_type => 1,
        },
        use_best_n_alleles => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 7,
            store       => \$use_best_n_alleles,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{freebayes};

    ## Options
    if ($apply_standard_filter) {

        # Equivalent to -m 30 -q 20 -R 0 -S 0
        push @commands, q{--standard-filters};
    }

    if ($calculate_genotype_quality) {

        push @commands, q{--genotype-qualities};
    }

    push @commands, q{--use-best-n-alleles} . $SPACE . $use_best_n_alleles;

    # Reference sequence file
    push @commands, q{--fasta-reference} . $SPACE . $referencefile_path;

    ## Infile
    push @commands, join $SPACE, @{$infile_paths_ref};

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
            stdoutfile_path => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
