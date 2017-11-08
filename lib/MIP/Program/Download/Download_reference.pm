package MIP::Program::Download::Download_reference;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ canonpath };

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
    our @EXPORT_OK = qw{ download_reference };
}

## Constants
Readonly my $SPACE => q{ };

sub download_reference {

## Function : Perl wrapper for writing download_reference command to filehandle. 
##          : Based on download_reference version 0.0.3
## Returns  : @commands
## Arguments: $reference_genome_versions_ref => Array with genome versions to downlaod {REF}
##          : $reference_dir_path            => Reference directory
##          : $stdoutfile_path               => Stdoutfile path
##          : $stderrfile_path               => Stderrfile path
##          : $stderrfile_path_append        => Append stderr info to file path
##          : $FILEHANDLE                    => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_genome_versions_ref;
    my $reference_dir_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        reference_genome_versions_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$reference_genome_versions_ref,
        },
        reference_dir_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$reference_dir_path,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{download_reference};

    push @commands,
      q{--reference_dir} . $SPACE . canonpath($reference_dir_path);

    push @commands,
        q{--reference_genome_versions}
      . $SPACE
      . join $SPACE
      . q{--reference_genome_versions}
      . $SPACE,
      @{$reference_genome_versions_ref};

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
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

    return @commands;
}

1;
