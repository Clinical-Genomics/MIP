package MIP::QC::Record;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    # Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];

use FindBin qw($Bin);    # Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );

BEGIN {
    use base qw (Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(add_program_info_to_sample_info);

}

sub add_program_info_to_sample_info {

##add_program_info_to_sample_info

##Function : Adds path and/or outdirectory and/or outfile from programs to sample_info to track all files and extract downstream
##Returns  : ""
##Arguments: $sample_info_href, $program_name, $path, $outdirectory, $outfile, $sample_id, $infile,
##         : $sample_info_href => Records on samples and family hash {REF}
##         : $program_name     => Program name
##         : $path             => Path of file
##         : $outdirectory     => Outdirectory of the file
##         : $outfile          => Outfile name
##         : $sample_id        => Sample_id for data at sample level {Optional}
##         : $infile           => Infile for data at sample level {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $program_name;
    my $path;
    my $outdirectory;
    my $outfile;
    my $sample_id;
    my $infile;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        path => {
            strict_type => 1,
            store       => \$path
        },
        outdirectory => {
            strict_type => 1,
            store       => \$outdirectory
        },
        outfile => {
            strict_type => 1,
            store       => \$outfile
        },
        infile    => { strict_type => 1, store => \$infile },
        sample_id => { strict_type => 1, store => \$sample_id },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    ## Sample level
    if ( defined $sample_id && defined $infile ) {

        if ( defined $path ) {

            # Path of file and/or directory
            $sample_info_href->{sample}{$sample_id}{program}{$program_name}
              {$infile}{path} = $path;
        }
        if ( defined $outdirectory ) {

            # Out directory of file
            $sample_info_href->{sample}{$sample_id}{program}{$program_name}
              {$infile}{outdirectory} = $outdirectory;
        }
        if ( defined $outfile ) {

            # Outfile
            $sample_info_href->{sample}{$sample_id}{program}{$program_name}
              {$infile}{outfile} = $outfile;

        }
    }
    else {

        ## Family level info
        if ( defined $path ) {

            $sample_info_href->{program}{$program_name}{path} = $path;
        }
        if ( defined $outdirectory ) {

            # Out directory of file
            $sample_info_href->{program}{$program_name}{outdirectory} =
              $outdirectory;
        }
        if ( defined $outfile ) {

            # Outfile
            $sample_info_href->{program}{$program_name}{outfile} = $outfile;
        }
    }
    return;
}

1;
