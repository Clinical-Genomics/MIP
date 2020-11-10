package MIP::Program::Upd;

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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ upd_call };
}

sub upd_call {

## Function : Perl wrapper for upd calls. Based on upd version 0.1
## Returns  : @commands
## Arguments: $af_tag                 => Allele tag to use
##          : $call_type              => Output regions or sites
##          : $father_id              => Father ID in vcf
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Path to infile vcf
##          : $mother_id              => Mother ID in vcf
##          : $outfile_path           => Path to outfile
##          : $proband_id             => Proband ID in vcf
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $af_tag;
    my $call_type;
    my $father_id;
    my $filehandle;
    my $infile_path;
    my $mother_id;
    my $outfile_path;
    my $proband_id;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    my $tmpl = {
        af_tag => {
            store       => \$af_tag,
            strict_type => 1,
        },
        call_type => {
            allow       => [qw{ regions sites }],
            defined     => 1,
            required    => 1,
            store       => \$call_type,
            strict_type => 1,
        },
        father_id => {
            defined     => 1,
            required    => 1,
            store       => \$father_id,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        mother_id => {
            defined     => 1,
            required    => 1,
            store       => \$mother_id,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        proband_id => {
            defined     => 1,
            required    => 1,
            store       => \$proband_id,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => q{upd}, } ), );

    if ($af_tag) {
        push @commands, q{--af-tag} . $SPACE . $af_tag;
    }

    push @commands, q{--vcf} . $SPACE . $infile_path;

    push @commands, q{--proband} . $SPACE . $proband_id;

    push @commands, q{--mother} . $SPACE . $mother_id;

    push @commands, q{--father} . $SPACE . $father_id;

    push @commands, $call_type;

    if ($outfile_path) {
        push @commands, q{--out} . $SPACE . $outfile_path;
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
