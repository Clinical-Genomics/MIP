package MIP::Program::Alignment::Chanjo;

use 5.026;
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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ chanjo_sex };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub chanjo_sex {

## Function : Perl wrapper for writing chanjo sex recipe to $filehandle. Based on chanjo 4.0.0
## Returns  : @commands
## Arguments: $chr_prefix                           => Chromosome prefix
##          : $filehandle                           => Sbatch filehandle to write to
##          : $infile_path                          => Infile path
##          : $log_file_path                        => Log file path
##          : $log_level                            => Level of logging
##          : $outfile_path                         => Outfile path
##          : $stderrfile_path                      => Stderrfile path
##          : $stderrfile_path_append               => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chr_prefix;
    my $filehandle;
    my $infile_path;
    my $log_file_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $log_level;

    my $tmpl = {
        chr_prefix => {
            allow       => [ undef, qw{chr} ],
            strict_type => 1,
            store       => \$chr_prefix
        },
        filehandle  => { required => 1, store => \$filehandle },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        log_file_path => { strict_type => 1, store => \$log_file_path },
        log_level     => {
            default     => q{INFO},
            allow       => [qw{DEBUG INFO WARNING ERROR CRITICAL}],
            strict_type => 1,
            store       => \$log_level
        },
        outfile_path           => { strict_type => 1, store => \$outfile_path },
        stderrfile_path        => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append => { strict_type => 1, store => \$stderrfile_path_append },
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

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            filehandle   => $filehandle,
        }
    );

    return @commands;
}

1;
