package MIP::Program::Gens;

use 5.026;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $ASTERISK $AMPERSAND $COLON $DOT $DOUBLE_QUOTE $EMPTY_STR $ESCAPE $NEWLINE $SPACE $UNDERSCORE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gens_generatedata };
}

sub gens_generatedata {

## Function : Perl wrapper for writing generate_gens_data.pl recipe to $filehandle.
## Returns  : @commands
## Arguments: $filehandle                                    => Sbatch filehandle to write to
##          : $gnomad_positions_ref                          => 
##          : $infile_tsv_path                               => Infile path
##          : $infile_vcf_path                               => Infile path
##          : $outfile_prefix                                => Outfile prefix
##          : $stderrfile_path                               => Stderrfile path
##          : $verbosity                                     => Set the minimum level of logging
##          : $xargs_mode                                    => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $gnomad_positions_ref;
    my $infile_tsv_path;
    my $infile_vcf_path;
    my $outfile_prefix;
    my $stderrfile_path;

    ## Default(s)
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        gnomad_positions_ref => {
            defined     => 1,
            required    => 1,
            store       => \$gnomad_positions_ref,
            strict_type => 1,
        },
        infile_tsv_path => {
            allow       => qr/ (?: tsv )$ /xms,
            defined     => 1,
            required    => 1,
            store       => \$infile_tsv_path,
            strict_type => 1,
        },
        infile_vcf_path => {
            allow       => qr/ (?: vcf.gz )$ /xms,
            defined     => 1,
            required    => 1,
            store       => \$infile_vcf_path,
            strict_type => 1,
        },
        outfile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_prefix,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => q{generate_gens_data.pl}, } ), );

    push @commands, $infile_tsv_path;
    push @commands, $infile_vcf_path;
    push @commands, $outfile_prefix;
    push @commands, $gnomad_positions_ref;

    push @commands,
    unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
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
