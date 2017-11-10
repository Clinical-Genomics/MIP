package MIP::Program::Variantcalling::Plink;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ plink plink_variant_pruning };
}

## Constants
Readonly my $SPACE => q{ };
Readonly my $COMMA => q{,};

sub plink_variant_pruning {

## Function : Perl wrapper for writing Plink recipe to prune variants and create unique IDs. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $vcffile_path        => Vcf file path
##          : $outfile_prefix      => Outfile prefix
##          : $stderrfile_path     => Stderrfile path
##          : $vcf_require_gt      => Skip variants where the GT field is absent
##          : $vcf_half_call       => Specify how '0/.' and similar VCF GT values should be handled
##          : $set_missing_var_ids => Assign chromosome-and-position-based IDs
##          : $const_fid           => Assign family id to all samples in file
##          : $make_bed            => Save data in plink binary format
##          : $indep               => Produce a pruned subset of markers that are in approximate linkage equilibrium with each other
##          : $indep_window_size   => Indep window size (kb)
##          : $indep_step_size     => Indep step size (variant ct)
##          : $indep_vif_threshold => Indep vif threshold
##          : $FILEHANDLE          => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $outfile_prefix;
    my $stderrfile_path;
    my $vcffile_path;
    my $vcf_half_call;
    my $set_missing_var_ids;
    my $const_fid;
    my $make_bed;
    my $indep;
    my $indep_window_size;
    my $indep_step_size;
    my $indep_vif_threshold;
    my $FILEHANDLE;

    ## Default(s)
    my $vcf_require_gt;

    my $tmpl = {
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        vcffile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$vcffile_path,
        },
        outfile_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_prefix,
        },
        vcf_require_gt => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$vcf_require_gt
        },
        vcf_half_call => {
            allow       => [ undef, qw{ error haploid missing reference } ],
            strict_type => 1,
            store       => \$vcf_half_call
        },
        set_missing_var_ids => {
            required    => 1,
            strict_type => 1,
            store       => \$set_missing_var_ids
        },
        const_fid => {
          required    => 1,
          strict_type => 1, store => \$const_fid },
        make_bed  => {
            strict_type => 1,
            store       => \$make_bed
        },
        indep => {
            required    => 1,
            strict_type => 1,
            store       => \$indep
        },
        indep_window_size =>
          { required => 1, strict_type => 1, store => \$indep_window_size },
        indep_step_size =>
          { required => 1, strict_type => 1, store => \$indep_step_size },
        indep_vif_threshold =>
          { required => 1, strict_type => 1, store => \$indep_vif_threshold },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{plink2};

    push @commands, q{--vcf} . $SPACE . $vcffile_path;

    if ($vcf_require_gt) {

        push @commands, q{--vcf-require-gt};
    }

    if ($vcf_half_call) {

        push @commands, q{--vcf-half-call} . $SPACE . $vcf_half_call;
    }

    push @commands, q{--set-missing-var-ids} . $SPACE . $set_missing_var_ids;

    push @commands, q{--const-fid} . $SPACE . $const_fid;

    if ($make_bed) {

        push @commands, q{--make-bed};
    }

    push @commands, q{--indep};
    push @commands, $indep_window_size;
    push @commands, $indep_step_size;
    push @commands, $indep_vif_threshold;


    ## Outfile
    push @commands, q{--out} . $SPACE . $outfile_prefix;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return;
}









1;
