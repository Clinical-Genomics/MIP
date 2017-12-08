package MIP::Program::Alignment::Chanjo;

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
use MIP::Script::Setup_script
  qw{ write_return_to_conda_environment write_source_environment_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ chanjo_sex };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub chanjo_sex {

## Function : Perl wrapper for writing chanjo sex recipe to $FILEHANDLE. Based on chanjo 4.0.0
## Returns  : @commands
## Arguments: $chr_prefix                           => Chromosome prefix
##          : $deactive_program_source              => Deactivate program specific environment
##          : $FILEHANDLE                           => Sbatch filehandle to write to
##          : $infile_path                          => Infile path
##          : $log_file_path                        => Log file path
##          : $log_level                            => Level of logging
##          : $outfile_path                         => Outfile path
##          : $program_source_command               => Program specific source environment commmand
##          : $source_main_environment_commands_ref => Source main environment command {REF}
##          : $stderrfile_path                      => Stderrfile path
##          : $stderrfile_path_append               => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chr_prefix;
    my $deactive_program_source;
    my $FILEHANDLE;
    my $infile_path;
    my $log_file_path;
    my $outfile_path;
    my $program_source_command;
    my $source_main_environment_commands_ref;
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
        deactive_program_source => {
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$deactive_program_source
        },
        FILEHANDLE  => { required => 1, store => \$FILEHANDLE },
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
        outfile_path => { strict_type => 1, store => \$outfile_path },
        program_source_command =>
          { strict_type => 1, store => \$program_source_command },
        source_main_environment_commands_ref => {
            default     => [],
            strict_type => 1,
            store       => \$source_main_environment_commands_ref,
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Write source program specific environment
    if ($program_source_command) {

        write_source_environment_command(
            {
                FILEHANDLE                      => $FILEHANDLE,
                source_environment_commands_ref => [$program_source_command],
            }
        );
    }

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
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    if ($deactive_program_source) {

        say {$FILEHANDLE} $NEWLINE;
        write_return_to_conda_environment(
            {
                FILEHANDLE => $FILEHANDLE,
                source_main_environment_commands_ref =>
                  $source_main_environment_commands_ref,
            }
        );
    }

    return @commands;
}

1;
