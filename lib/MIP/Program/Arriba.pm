package MIP::Program::Arriba;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $EQUALS $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    our @EXPORT_OK = qw{ arriba draw_fusions };
}

sub arriba {

## Function : Perl wrapper for Arriba commands module. Based on Arriba 1.1.0
## Returns  : @commands
## Arguments: $annotation_file_path       => Path to GTF file with annotations
##          : $blacklist_file_path        => Path to file with blacklist events
##          : $dna_sv_file_path           => Path to file with wgs SV calls
##          : $discarded_fusion_file_path => Path to write discarded fusion events to
##          : $filehandle                 => Filehandle to write to
##          : $genome_file_path           => Genome reference path
##          : $infile_path                => SAM/BAM file path
##          : $known_fusion_file_path     => Path to file with known fusions
##          : $outfile_path               => Path to outfile
##          : $protein_domain_file_path   => Path to file with protein domains
##          : $stderrfile_path            => Stderrfile path
##          : $stderrfile_path_append     => Append stderr info to file path
##          : $stdinfile_path             => Stdinfile path
##          : $stdoutfile_path            => Stdoutfile path
##          : $tag_file_path              => Tag file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotation_file_path;
    my $blacklist_file_path;
    my $dna_sv_file_path;
    my $discarded_fusion_file_path;
    my $filehandle;
    my $genome_file_path;
    my $infile_path;
    my $known_fusion_file_path;
    my $outfile_path;
    my $protein_domain_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $tag_file_path;

    my $tmpl = {
        annotation_file_path => {
            required    => 1,
            store       => \$annotation_file_path,
            strict_type => 1,
        },
        blacklist_file_path => {
            store       => \$blacklist_file_path,
            strict_type => 1,
        },
        dna_sv_file_path => {
            store       => \$dna_sv_file_path,
            strict_type => 1,
        },
        discarded_fusion_file_path => {
            store       => \$discarded_fusion_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        genome_file_path => {
            required    => 1,
            store       => \$genome_file_path,
            strict_type => 1,
        },
        infile_path => {
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        known_fusion_file_path => {
            store       => \$known_fusion_file_path,
            strict_type => 1,
        },
        outfile_path => {
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        protein_domain_file_path => {
            store       => \$protein_domain_file_path,
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
        stdinfile_path => {
            store       => \$stdinfile_path,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        tag_file_path => {
            store       => \$tag_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => q{arriba}, } ), );

    push @commands, q{-g} . $SPACE . $annotation_file_path;

    if ($blacklist_file_path) {

        push @commands, q{-b} . $SPACE . $blacklist_file_path;
    }

    if ($dna_sv_file_path) {

        push @commands, q{-d} . $SPACE . $dna_sv_file_path;
    }

    if ($discarded_fusion_file_path) {

        push @commands, q{-O} . $SPACE . $discarded_fusion_file_path;
    }

    push @commands, q{-a} . $SPACE . $genome_file_path;

    push @commands, q{-x} . $SPACE . $infile_path;

    if ($known_fusion_file_path) {

        push @commands, q{-k} . $SPACE . $known_fusion_file_path;
    }

    push @commands, q{-o} . $SPACE . $outfile_path;

    if ($protein_domain_file_path) {

        push @commands, q{-p} . $SPACE . $protein_domain_file_path;
    }

    if ($tag_file_path) {

        push @commands, q{-t} . $SPACE . $tag_file_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub draw_fusions {

## Function : Perl wrapper for Arriba's draw_fusion.R script. Based on Arriba 1.1.0
## Returns  : @commands
## Arguments: $alignment_file_path      => Path to BAM file with alignments
##          : $annotation_file_path     => Path to annotation file
##          : $cytoband_file_path       => Path to file with coordinates for cytobands
##          : $filehandle               => Filehandle to write to
##          : $fusion_file_path         => Path gene fusion file
##          : $outfile_path             => Path to PDF outfile
##          : $protein_domain_file_path => Path to file with protein domains
##          : $stderrfile_path          => Stderrfile path
##          : $stderrfile_path_append   => Append stderr info to file path
##          : $stdoutfile_path          => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $alignment_file_path;
    my $annotation_file_path;
    my $cytoband_file_path;
    my $filehandle;
    my $fusion_file_path;
    my $outfile_path;
    my $protein_domain_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        alignment_file_path => {
            required    => 1,
            store       => \$alignment_file_path,
            strict_type => 1,
        },
        annotation_file_path => {
            required    => 1,
            store       => \$annotation_file_path,
            strict_type => 1,
        },
        cytoband_file_path => {
            store       => \$cytoband_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        fusion_file_path => {
            required    => 1,
            store       => \$fusion_file_path,
            strict_type => 1,
        },
        outfile_path => {
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        protein_domain_file_path => {
            store       => \$protein_domain_file_path,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => q{draw_fusions.R}, } ), );

    push @commands, q{--alignments} . $EQUALS . $alignment_file_path;

    push @commands, q{--annotation} . $EQUALS . $annotation_file_path;

    if ($cytoband_file_path) {

        push @commands, q{--cytobands} . $EQUALS . $cytoband_file_path;
    }

    push @commands, q{--fusions} . $EQUALS . $fusion_file_path;

    push @commands, q{--output} . $EQUALS . $outfile_path;

    if ($protein_domain_file_path) {

        push @commands, q{--proteinDomains} . $EQUALS . $protein_domain_file_path;
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
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}
1;
