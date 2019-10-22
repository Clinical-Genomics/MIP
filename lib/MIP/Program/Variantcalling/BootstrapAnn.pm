package MIP::Program::Variantcalling::BootstrapAnn;

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
    our @EXPORT_OK = qw{ bootstrapann };
}

## Constants
Readonly my $SPACE => q{ };

sub bootstrapann {

## Function : Perl wrapper for generic commands module.
## Returns  : @commands
## Arguments: $ase_file_path          => Path to file with GATK ASEReadCounter data
##          : $filehandle             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $vcf_infile_path        => VCF infile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $ase_file_path;
    my $filehandle;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $vcf_infile_path;

    ## Default(s)

    my $tmpl = {
        ase_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$ase_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
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
        vcf_infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$vcf_infile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{BootstrapAnn.py};

    push @commands, q{--vcf} . $SPACE . $vcf_infile_path;

    push @commands, q{--ase} . $SPACE . $ase_file_path;

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
            filehandle   => $filehandle,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
