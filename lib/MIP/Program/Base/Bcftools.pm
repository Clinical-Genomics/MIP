package MIP::Program::Base::Bcftools;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Params::Check qw{ check allow last_error };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ bcftools_base };
}

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

sub bcftools_base {

## Function : Perl wrapper for picardtools base. Based on Picardtools v2.9.2-SNAPSHOT
## Returns  : @commands
## Arguments: $commands_ref => List of commands added earlier
##          : $FILEHANDLE   => Filehandle to write to
##          : $outfile_path => Outfile path
##          : $output_type  => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref  => Regions to process {REF}
##          : $samples_file => File of samples to annotate

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $FILEHANDLE;
    my $outfile_path;
    my $regions_ref;
    my $samples_file;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        commands_ref =>
          { default => [], strict_type => 1, store => \$commands_ref },
        FILEHANDLE   => { store       => \$FILEHANDLE },
        outfile_path => { strict_type => 1, store => \$outfile_path, },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v}],
            strict_type => 1,
            store       => \$output_type,
        },
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref, },
        samples_file => { strict_type => 1, store => \$samples_file, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = @{$commands_ref};

    if ($samples_file) {

        push @commands, q{--samples-file} . $SPACE . $samples_file;
    }
    if ( @{$regions_ref} ) {

        # Limit output to regions
        push @commands, q{--regions_ref} . $SPACE . join $COMMA,
          @{$regions_ref};
    }
    if ($outfile_path) {

        # Specify output filename
        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    #Specify output type
    push @commands, q{--output-type} . $SPACE . $output_type;

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
