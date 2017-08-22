package MIP::Program::Alignment::Chanjo;

use Carp;
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };

use FindBin qw{$Bin};    #Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    use base qw {Exporter};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{chanjo_sex};

}

use Params::Check qw{check allow last_error};
use Readonly;

Readonly my $SPACE => q{ };

sub chanjo_sex {

##chanjo_sex

##Function : Perl wrapper for writing chanjo sex recipe to $FILEHANDLE. Based on chanjo 4.0.0
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $log_file_path, $chr_prefix, $log_level
##         : $FILEHANDLE      => Sbatch filehandle to write to
##         : $infile_path     => Infile path
##         : $outfile_path    => Outfile path
##         : $stderrfile_path => Stderrfile path
##         : $log_file_path   => Log file path
##         : $chr_prefix      => Chromosome prefix
##         : $log_level       => Level of logging

    my ($arg_href) = @_;

    ## Default(s)
    my $log_level;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $log_file_path;
    my $chr_prefix;
    my $stderrfile_path_append;

    my $tmpl = {
        FILEHANDLE  => { required => 1, store => \$FILEHANDLE },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        log_file_path   => { strict_type => 1, store => \$log_file_path },
        chr_prefix      => {
            allow       => [qw{undef chr}],
            strict_type => 1,
            store       => \$chr_prefix
        },
        log_level => {
            default     => q{INFO},
            allow       => [qw{DEBUG INFO WARNING ERROR CRITICAL}],
            strict_type => 1,
            store       => \$log_level
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Chanjo
    my @commands = qw{chanjo};    #Stores commands depending on input parameters

    ## Chanjo main options
    if ($log_level) {

        push @commands, q{--log-level} . $SPACE . $log_level;
    }
    if ($log_file_path) {

        push @commands, q{--log-file} . $SPACE . $log_file_path;
    }

    push @commands, q{sex};

    ## Options
    if ($chr_prefix) {

        push @commands, q{--prefix} . $SPACE . $chr_prefix;
    }
    ##Infile
    push @commands, $infile_path;

    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }
    if ($stderrfile_path) {

        push @commands,
          unix_standard_streams(
            {
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
            }
          );
    }
    if ($FILEHANDLE) {

        unix_write_to_file(
            {
                commands_ref => \@commands,
                separator    => $SPACE,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    return @commands;
}

1;
