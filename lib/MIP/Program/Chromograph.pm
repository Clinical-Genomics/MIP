package MIP::Program::Chromograph;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ chromograph };
}

sub chromograph {

## Function : Perl wrapper for chromograph.
## Returns  : @commands
## Arguments: $autozyg_file_path      => Autozygosity data, bed file
##          : $coverage_file_path     => Coverage data infile
##          : $euploid                => Generate png files for all chromosomes
##          : $filehandle             => Filehandle to write to
##          : $fracsnp_file_path      => Fraction of homozygous SNP, wig file
##          : $ideogram_file_path     => Bed file with ideogram data
##          : $outdir_path            => Outdir path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $step                   => Bin size for wig file
##          : $upd_regions_file_path  => UPD regions data file
##          : $upd_sites_file_path    => UPD sites data file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $autozyg_file_path;
    my $coverage_file_path;
    my $euploid;
    my $filehandle;
    my $fracsnp_file_path;
    my $ideogram_file_path;
    my $outdir_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $step;
    my $upd_regions_file_path;
    my $upd_sites_file_path;

    my $tmpl = {
        autozyg_file_path => {
            allow       => qr/ [.] bed \z /xms,
            store       => \$autozyg_file_path,
            strict_type => 1,
        },
        coverage_file_path => {
            allow       => qr/ [.] wig \z /xms,
            store       => \$coverage_file_path,
            strict_type => 1,
        },
        euploid => {
            allow       => [ undef, 0, 1 ],
            store       => \$euploid,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        fracsnp_file_path => {
            allow       => qr/ [.] wig \z /xms,
            store       => \$fracsnp_file_path,
            strict_type => 1,
        },
        ideogram_file_path => {
            allow       => qr/ [.] bed \z /xms,
            store       => \$ideogram_file_path,
            strict_type => 1,
        },
        outdir_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path,
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
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        step => {
            allow       => qr/\A \d+ \z/xms,
            store       => \$step,
            strict_type => 1,
        },
        upd_regions_file_path => {
            allow       => qr/ [.] bed \z /xms,
            store       => \$upd_regions_file_path,
            strict_type => 1,
        },
        upd_sites_file_path => {
            allow       => qr/ [.] bed \z /xms,
            store       => \$upd_sites_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => q{chromograph}, } ), );

    if ($autozyg_file_path) {

        push @commands, q{--autozyg} . $SPACE . $autozyg_file_path;
    }

    if ($coverage_file_path) {

        push @commands, q{--coverage} . $SPACE . $coverage_file_path;
    }

    if ($euploid) {

        push @commands, q{--euploid};
    }

    if ($fracsnp_file_path) {

        push @commands, q{--fracsnp} . $SPACE . $fracsnp_file_path;
    }

    if ($ideogram_file_path) {

        push @commands, q{--ideogram} . $SPACE . $ideogram_file_path;
    }

    push @commands, q{--outd} . $SPACE . $outdir_path;

    if ($step) {

        push @commands, q{--step} . $SPACE . $step;
    }

    if ($upd_regions_file_path) {

        push @commands, q{--regions} . $SPACE . $upd_regions_file_path;
    }

    if ($upd_sites_file_path) {

        push @commands, q{--sites} . $SPACE . $upd_sites_file_path;
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

1;
