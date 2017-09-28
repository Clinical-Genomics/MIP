package MIP::Program::Variantcalling::Tiddit;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

use Readonly;

use FindBin qw{ $Bin };    #Find directory of script
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ tiddit_sv };

}

## Constants
Readonly my $SPACE => q{ };

sub tiddit_sv {

## tiddit_sv

## Function : Perl wrapper for writing tiddit sv recipe to $FILEHANDLE or return commands array. Based on tiddit 1.0.2.
## Returns  : "@commands"
## Arguments: $infile_path, $outfile_path_prefix, $FILEHANDLE, $minimum_number_supporting_pairs
##          : $infile_path                     => Infile path
##          : $outfile_path_prefix             => Outfile path. Write documents to FILE
##          : $FILEHANDLE                      => Filehandle to write to
##          : $minimum_number_supporting_pairs => Minimum number of supporting pairs in order to call a variation event

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path_prefix;
    my $FILEHANDLE;
    my $minimum_number_supporting_pairs;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path_prefix =>
          { strict_type => 1, store => \$outfile_path_prefix },
        FILEHANDLE                      => { store => \$FILEHANDLE },
        minimum_number_supporting_pairs => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$minimum_number_supporting_pairs
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## tiddit
    my @commands = qw{ TIDDIT --sv };

    ## Option: minimum number of supporting pairs in order to call a variation event
    if ($minimum_number_supporting_pairs) {

        push @commands, q{-p} . $SPACE . $minimum_number_supporting_pairs;
    }

    #Outfile prefix
    if ($outfile_path_prefix) {

        push @commands, q{-o} . $SPACE . $outfile_path_prefix;
    }

    ## Infile
    push @commands, q{-b} . $SPACE . $infile_path;

    push @commands,
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
