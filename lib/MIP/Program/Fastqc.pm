package MIP::Program::Fastqc;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ fastqc };

}

Readonly my $BASE_COMMAND => q{fastqc};

sub fastqc {

##Function : Perl wrapper for writing fastqc recipe to already open $filehandle or return commands array. Based on fastq 0.11.5
##Returns  : @commands
##Arguments: $extract           => If set then the zipped output file will be uncompressed in the same directory after it has been created
##         : $filehandle        => Filehandle to write to
##         : $infile_path       => Infile path
##         : $outdirectory_path => Outdirectory path
##         : $quiet             => Supress all progress messages on stdout and only report errors

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outdirectory_path;

    ## Default(s)
    my $extract;
    my $quiet;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outdirectory_path => { store => \$outdirectory_path, strict_type => 1, },
        extract           => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$extract,
            strict_type => 1,
        },
        quiet => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$quiet,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), );

    if ($quiet) {

        push @commands, q{--quiet};
    }
    if ($extract) {

        push @commands, q{--extract};
    }
    if ($outdirectory_path) {

        push @commands, q{--outdir} . $SPACE . $outdirectory_path;
    }

    push @commands, $infile_path;

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
