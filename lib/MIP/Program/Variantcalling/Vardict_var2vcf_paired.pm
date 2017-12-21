package MIP::Program::Variantcalling::Vardict_var2vcf_paired;

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
    our @EXPORT_OK = qw{ vardict_var2vcf_paired };
}

## Constants
Readonly my $DQUOTE => q{"};
Readonly my $PIPE   => q{|};
Readonly my $SPACE  => q{ };

sub vardict_var2vcf_paired {

## Function : Perl wrapper for var2vcf_paired, converting variant output to vcf. Based on vardict 2017.09.24
## Returns  : @commands
## Arguments: $af_threshold           => Threshold for allele frequency
##          : $FILEHANDLE             => Filehandle to write to
##          : $sample_name            => Sample name to be used directly, will overwrite -n option
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $sample_name;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $af_threshold;

    my $tmpl = {
        af_threshold => {
            ## Exactly 2 decimal points after 0 or 1
            allow       => qr/ ^0.\d{1,2}$ | ^1$ /xsm,
            defined  => 1,
            default  => 0.01,
            required => 1,
            store       => \$af_threshold,
            strict_type => 1,
        },
        sample_name => {
            allow       => qr/ ^\w+$ /xsm,
            required    => 1,
            store       => \$sample_name,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
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

    ## Stores commands depending on input parameters
    my @commands = q{var2vcf_paired.pl};

    push @commands, q{-f} . $SPACE . $af_threshold;

    push @commands, q{-N} . $SPACE . $sample_name;

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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
