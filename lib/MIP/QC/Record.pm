package MIP::QC::Record;

use strict;
use warnings;
use warnings qw{FATAL utf8};
use utf8;    # Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{check allow last_error};

use FindBin qw{$Bin};    # Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{add_program_outfile_to_sample_info add_program_metafile_to_sample_info add_processing_metafile_to_sample_info};

}

sub add_program_outfile_to_sample_info {

##add_program_outfile_to_sample_info

##Function : Adds path and/or outdirectory and/or outfile and/or version from programs to sample_info to track all outfiles and extract downstream
##Returns  : ""
##Arguments: $sample_info_href, $program_name, $path, $outdirectory, $outfile, $sample_id, $infile,
##         : $sample_info_href => Records on samples and family hash {REF}
##         : $program_name     => Program name
##         : $path             => Path of file
##         : $outdirectory     => Outdirectory of the file
##         : $outfile          => Outfile name
##         : $version          => Version of file
##         : $sample_id        => Sample_id for data at sample level {Optional}
##         : $infile           => Infile for data at sample level {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $program_name;
    my $path;
    my $outdirectory;
    my $outfile;
    my $version;
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
        version => {
            strict_type => 1,
            store       => \$version
        },
        sample_id => { strict_type => 1, store => \$sample_id },
        infile    => { strict_type => 1, store => \$infile },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set the key and value pair to add to sample_info hash
    my %parameter = (
        path         => $path,
        outdirectory => $outdirectory,
        outfile      => $outfile,
        version      => $version,
    );

    if ( defined $sample_id && defined $infile ) {

      SAMPLE_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{sample}{$sample_id}{program}{$program_name}
                  {$infile}{$parameter_key} = $parameter_value;
            }
        }
    }
    else {

      FAMILY_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{program}{$program_name}{$parameter_key} =
                  $parameter_value;
            }
        }
    }
    return;
}

sub add_program_metafile_to_sample_info {

##add_program_metafile_to_sample_info

##Function : Adds path and/or directory and/or file and/or version from programs to sample_info to track all metafiles and extract downstream
##Returns  : ""
##Arguments: $sample_info_href, $program_name, $metafile_tag, $path, $directory, $file, $version, $sample_id, $infile,
##         : $sample_info_href => Records on samples and family hash {REF}
##         : $program_name     => Program name
##         : $metafile_tag     => Id tag of meta file
##         : $path             => Path of file
##         : $directory        => Directory of the file
##         : $file             => File name
##         : $version          => Version of file
##         : $sample_id        => Sample_id for data at sample level {Optional}
##         : $infile           => Infile for data at sample level {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $program_name;
    my $metafile_tag;
    my $path;
    my $directory;
    my $file;
    my $version;
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
        metafile_tag => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$metafile_tag
        },
        path => {
            strict_type => 1,
            store       => \$path
        },
        directory => {
            strict_type => 1,
            store       => \$directory
        },
        file => {
            strict_type => 1,
            store       => \$file
        },
        version => {
            strict_type => 1,
            store       => \$version
        },
        sample_id => { strict_type => 1, store => \$sample_id },
        infile    => { strict_type => 1, store => \$infile },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set the key and value pair to add to sample_info hash
    my %parameter = (
        path      => $path,
        directory => $directory,
        file      => $file,
        version   => $version,
    );

    if ( defined $sample_id && defined $infile ) {

      SAMPLE_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{sample}{$sample_id}{program}{$program_name}
                  {$infile}{$metafile_tag}{$parameter_key} = $parameter_value;
            }
        }
    }
    else {

      FAMILY_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{program}{$program_name}{$metafile_tag}
                  {$parameter_key} = $parameter_value;
            }
        }
    }
    return;
}

sub add_processing_metafile_to_sample_info {

##add_processing_metafile_to_sample_info

##Function : Adds metafile path from sample_id|family_id processing to sample_info to track all metafiles and extract downstream
##Returns  : ""
##Arguments: $sample_info_href, $metafile_tag, $path, $sample_id
##         : $sample_info_href => Records on samples and family hash {REF}
##         : $metafile_tag     => Id tag of meta file
##         : $path             => Path of file
##         : $sample_id        => Sample_id for data at sample level {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $metafile_tag;
    my $path;
    my $sample_id;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        metafile_tag => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$metafile_tag
        },
        path => {
            strict_type => 1,
            store       => \$path
        },
        sample_id => { strict_type => 1, store => \$sample_id },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set the key and value pair to add to sample_info hash
    my %parameter = ( path => $path, );

    if ( defined $sample_id ) {

      SAMPLE_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{sample}{$sample_id}{$metafile_tag}
                  {$parameter_key} = $parameter_value;
            }
        }
    }
    else {

      FAMILY_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{$metafile_tag}{$parameter_key} =
                  $parameter_value;
            }
        }
    }
    return;
}

1;
