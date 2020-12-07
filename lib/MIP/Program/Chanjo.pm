package MIP::Program::Chanjo;

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
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ chanjo_sex };
}

Readonly my $BASE_COMMAND => q{chanjo};

sub chanjo_sex {

## Function : Perl wrapper for writing chanjo sex recipe to $filehandle. Based on chanjo 4.0.0
## Returns  : @commands
## Arguments: $chr_prefix             => Chromosome prefix
##          : $filehandle             => Sbatch filehandle to write to
##          : $infile_path            => Infile path
##          : $log_file_path          => Log file path
##          : $log_level              => Level of logging
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append

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
            allow       => [ undef, qw{ chr } ],
            store       => \$chr_prefix,
            strict_type => 1,
        },
        filehandle  => { required => 1, store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        log_file_path => { store => \$log_file_path, strict_type => 1, },
        log_level     => {
            allow       => [qw{ DEBUG INFO WARNING ERROR CRITICAL }],
            default     => q{INFO},
            store       => \$log_level,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path,    strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), );

    if ($log_level) {

        push @commands, q{--log-level} . $SPACE . $log_level;
    }
    if ($log_file_path) {

        push @commands, q{--log-file} . $SPACE . $log_file_path;
    }

    push @commands, q{sex};

    if ($chr_prefix) {

        push @commands, q{--prefix} . $SPACE . $chr_prefix;
    }

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
