package MIP::Program::Qc::Theta;

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
    our @EXPORT_OK = qw{ run_theta  };
}

## Constants
Readonly my $SPACE => q{ };

sub run_theta {

    ## Function : Perl wrapper for THetA, a tool for estimating tumor purity and clonal/subclonal copy number aberrations based on THetA2 Jul 31, 2017.
    ## Returns  : @commands
    ## Arguments: $FILEHANDLE             => Filehandle to write to
    ##          : $normal_file            => File containing allelic counts for normal sample SNPs
    ##          : $num_processes          => The number of processes to be used
    ##          : $output_dir             => Directory where result file is written to
    ##          : $output_prefix          => A prefix for the output files
    ##          : $stderrfile_path        => Stderrfile path
    ##          : $stderrfile_path_append => Append stderr info to file path
    ##          : $stdoutfile_path        => Stdoutfile path
    ##          : $tumor_file             => File containing allelic counts for tumor sample SNPs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $normal_file;
    my $num_processes;
    my $output_dir;
    my $output_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $tumor_file;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        normal_file => {
            required    => 1,
            store       => \$normal_file,
            strict_type => 1,
        },
        num_processes => {
            allow       => qr/ ^\d+$ /xsm,
            store       => \$num_processes,
            strict_type => 1,
        },
        output_dir => {
            store       => \$output_dir,
            strict_type => 1,
        },
        output_prefix => {
            allow       => qr/ ^\w+$ /xsm,
            store       => \$output_prefix,
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
        tumor_file => {
            required    => 1,
            store       => \$tumor_file,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{RunTHetA};

    if ($output_dir) {

        #output directory
        push @commands, q{--DIR} . $SPACE . $output_dir;
    }

    if ($output_prefix) {

        #output prefix
        push @commands, q{--OUTPUT_PREFIX} . $SPACE . $output_prefix;
    }

    if ($num_processes) {

        #number of processes
        push @commands, q{--NUM_PROCESSES} . $SPACE . $num_processes;
    }

    push @commands, q{--NORMAL_FILE} . $SPACE . $normal_file;

    push @commands, q{--TUMOR_FILE} . $SPACE . $tumor_file;

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
