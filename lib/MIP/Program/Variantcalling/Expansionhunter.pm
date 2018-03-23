package MIP::Program::Variantcalling::Expansionhunter;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
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
    our @EXPORT_OK = qw{ expansionhunter };
}

## Constants
Readonly my $SPACE => q{ };

sub expansionhunter {

## Function : Perl wrapper for writing Expansion Hunter command to $FILEHANDLE or return commands array.
## Returns  : @commands
## Arguments: $FILEHANDLE              => Filehandle to write to
##          : $infile_path             => Path to sorted, indexed bam file
##          : $json_outfile_path       => Path to json file
##          : $log_outfile_path        => Path to log file
##          : $min_anchor_mapq         => Min MAPQ of an in-repeat read anchor
##          : $min_baseq               => Min base quality of a high confidece base call
##          : $min_score               => Min weigheted matching score
##          : $read_depth              => Expansion Hunter will not calculate read depth if set
##          : $reference_genome_path   => Path to reference fasta file
##          : $region_extension_length => Specifies how far from on/off-target regions to search for informative reads
##          : $repeat_specs_dir_path   => Path to dir with repeat-specification files
##          : $sex                     => Sex of the sample (Defaults to female)
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Append stderr info to file path
##          : $stdoutfile_path         => Stdoutfile path
##          : $vcf_outfile_path        => Path to vcf file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $json_outfile_path;
    my $log_outfile_path;
    my $min_anchor_mapq;
    my $min_baseq;
    my $min_score;
    my $read_depth;
    my $reference_genome_path;
    my $region_extension_length;
    my $repeat_specs_dir_path;
    my $sex;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $vcf_outfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        json_outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$json_outfile_path,
            strict_type => 1,
        },
        log_outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$log_outfile_path,
            strict_type => 1,
        },
        min_anchor_mapq => {
            allow       => qr/ ^\d+$ /xms,
            defined     => 1,
            store       => \$min_anchor_mapq,
            strict_type => 1,
        },
        min_baseq => {
            allow       => qr/ ^\d+$ /xms,
            defined     => 1,
            store       => \$min_baseq,
            strict_type => 1,
        },
        min_score => {
            allow       => qr/ ^0[.]\d+$ /xms,
            defined     => 1,
            store       => \$min_score,
            strict_type => 1,
        },
        read_depth => {
            allow       => qr/ ^\d+$ /xms,
            defined     => 1,
            store       => \$read_depth,
            strict_type => 1,
        },
        reference_genome_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_path,
            strict_type => 1,
        },
        region_extension_length => {
            allow       => qr/ ^\d+$ /xms,
            defined     => 1,
            store       => \$region_extension_length,
            strict_type => 1,
        },
        repeat_specs_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$repeat_specs_dir_path,
            strict_type => 1,
        },
        sex => {
            allow       => [qw{ female male }],
            defined     => 1,
            store       => \$sex,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        vcf_outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$vcf_outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{ExpansionHunter};

    ## Required arguments
    push @commands, q{--bam} . $SPACE . $infile_path;

    push @commands, q{--json} . $SPACE . $json_outfile_path;

    push @commands, q{--log} . $SPACE . $log_outfile_path;

    push @commands, q{--ref-fasta} . $SPACE . $reference_genome_path;

    push @commands, q{--repeat-specs} . $SPACE . $repeat_specs_dir_path;

    push @commands, q{--vcf} . $SPACE . $vcf_outfile_path;

    ## Optional arguments
    if ($min_anchor_mapq) {
        push @commands, q{--min-anchor-mapq} . $SPACE . $min_anchor_mapq;
    }

    if ($min_baseq) {
        push @commands, q{--min-baseq} . $SPACE . $min_baseq;
    }

    if ($min_score) {
        push @commands, q{--min-score} . $SPACE . $min_score;
    }

    if ($read_depth) {
        push @commands, q{--read-depth} . $SPACE . $read_depth;
    }

    if ($region_extension_length) {
        push @commands,
          q{--region-extension-length} . $SPACE . $region_extension_length;
    }

    if ($sex) {
        push @commands, q{--sex} . $SPACE . $sex;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
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
