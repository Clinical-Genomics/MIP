package MIP::Program::Vcf2cytosure;

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
    our @EXPORT_OK = qw{ vcf2cytosure_convert };
}

sub vcf2cytosure_convert {

## Function : Perl wrapper for Vcf2cytosure 0.5.0.
## Returns  : @commands
## Arguments: $blacklist_file_path    => List of regions to exclude file path
##          : $coverage_file          => Path to coverage file
##          : $filehandle             => Filehandle to write to
##          : $frequency              => Maximum frequency
##          : $frequency_tag          => Frequency tag of the info field
##          : $maxbnd                 => Maximum BND size
##          : $no_filter              => Disable any filtering
##          : $outfile_path           => Outfile path to write to
##          : $sex                    => Sex of sample
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $variant_size           => Minimum variant size.
##          : $vcf_infile_path        => Path to VCF infile
##          : $version                => Version of program

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $blacklist_file_path;
    my $coverage_file;
    my $filehandle;
    my $maxbnd;
    my $outfile_path;
    my $sex;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $vcf_infile_path;

    ## Default(s)
    my $frequency;
    my $frequency_tag;
    my $no_filter;
    my $variant_size;
    my $version;

    my $tmpl = {
        blacklist_file_path => {
            store       => \$blacklist_file_path,
            strict_type => 1,
        },
        coverage_file => {
            defined     => 1,
            required    => 1,
            store       => \$coverage_file,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        frequency => {
            default     => 0.01,
            store       => \$frequency,
            strict_type => 1,
        },
        frequency_tag => {
            default     => q{FRQ},
            store       => \$frequency_tag,
            strict_type => 1,
        },
        maxbnd => {
            allow       => qr/ ^\d+$ /xsm,
            store       => \$maxbnd,
            strict_type => 1,
        },
        no_filter => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$no_filter,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        sex => {
            allow       => [ undef, qw{ female male } ],
            store       => \$sex,
            strict_type => 1,
        },
        sex => {
            allow       => [ undef, qw{ female male } ],
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
        variant_size => {
            allow       => qr/ \A \d+ \z /xsm,
            default     => 5000,
            store       => \$variant_size,
            strict_type => 1,
        },
        vcf_infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$vcf_infile_path,
            strict_type => 1,
        },
        version => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => q{vcf2cytosure}, } ), );

    if ($version) {

        push @commands, q{--version};
    }

    if ($blacklist_file_path) {

        push @commands, q{--blacklist} . $SPACE . $blacklist_file_path;
    }
    if ($variant_size) {

        push @commands, q{--size} . $SPACE . $variant_size;
    }

    if ($frequency) {

        push @commands, q{--frequency} . $SPACE . $frequency;
    }

    if ($frequency_tag) {

        push @commands, q{--frequency_tag} . $SPACE . $frequency_tag;
    }

    if ($maxbnd) {

        push @commands, q{--maxbnd} . $SPACE . $maxbnd;
    }

    if ($no_filter) {

        push @commands, q{--no-filter};
    }

    if ($sex) {

        push @commands, q{--sex} . $SPACE . $sex;
    }

    if ($outfile_path) {

        push @commands, q{--out} . $SPACE . $outfile_path;
    }

    if ($sex) {

        push @commands, q{--sex} . $SPACE . $sex;
    }

    push @commands, q{--coverage} . $SPACE . $coverage_file;

    push @commands, q{--vcf} . $SPACE . $vcf_infile_path;

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
