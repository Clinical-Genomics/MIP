package MIP::QC::Record;

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

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ add_most_complete_vcf add_program_outfile_to_sample_info add_processing_metafile_to_sample_info add_program_metafile_to_sample_info };

}

sub add_most_complete_vcf {

## Function : Adds the most complete vcf file to sample_info
## Returns  :
## Arguments: $active_parameter_href     => Active parameters for this analysis hash {REF}
##          : $path                      => Path to file
##          : $program_name              => Program name
##          : $sample_info_href          => Info on samples and family hash {REF}
##          : $vcf_file_key              => Key for labelling most complete vcf
##          : $vcfparser_outfile_counter => Number of outfile files from in vcfParser (select, range)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $path;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $vcf_file_key;
    my $vcfparser_outfile_counter;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        path =>
          { required => 1, defined => 1, strict_type => 1, store => \$path },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        vcf_file_key => {
            default     => q{vcf_file},
            allow       => [qw{ vcf_file sv_vcf_file}],
            strict_type => 1,
            store       => \$vcf_file_key
        },
        vcfparser_outfile_counter => {
            default     => 0,
            strict_type => 1,
            store       => \$vcfparser_outfile_counter
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $active_parameter_href->{ q{p} . $program_name } == 1 ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            $sample_info_href->{$vcf_file_key}{clinical}{path} = $path;
        }
        else {

            $sample_info_href->{$vcf_file_key}{research}{path} = $path;
        }
    }
    return;
}

sub add_program_outfile_to_sample_info {

## Function : Adds path and/or outdirectory and/or outfile and/or version from programs to sample_info to track all outfiles and extract downstream
## Returns  :
## Arguments: $infile           => Infile for data at sample level {Optional}
##          : $outdirectory     => Outdirectory of the file
##          : $outfile          => Outfile name
##          : $path             => Path of file
##          : $program_name     => Program name
##          : $sample_id        => Sample_id for data at sample level {Optional}
##          : $sample_info_href => Records on samples and family hash {REF}
##          : $version          => Version of file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile;
    my $outdirectory;
    my $outfile;
    my $path;
    my $program_name;
    my $sample_id;
    my $sample_info_href;
    my $version;

    my $tmpl = {
        infile       => { strict_type => 1, store => \$infile },
        outdirectory => {
            strict_type => 1,
            store       => \$outdirectory,
        },
        outfile => {
            strict_type => 1,
            store       => \$outfile,
        },
        path => {
            strict_type => 1,
            store       => \$path,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_id        => { strict_type => 1, store => \$sample_id, },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        version => {
            strict_type => 1,
            store       => \$version,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set the key and value pair to add to sample_info hash
    my %parameter = (
        outdirectory => $outdirectory,
        outfile      => $outfile,
        path         => $path,
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

sub add_processing_metafile_to_sample_info {

## Function : Adds metafile path from sample_id|family_id processing to sample_info to track all metafiles and extract downstream
## Returns  :
## Arguments: $metafile_tag     => Id tag of meta file
##          : $path             => Path of file
##          : $sample_id        => Sample_id for data at sample level {Optional}
##          : $sample_info_href => Records on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $metafile_tag;
    my $path;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        metafile_tag => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$metafile_tag,
        },
        path => {
            strict_type => 1,
            store       => \$path,
        },
        sample_id        => { strict_type => 1, store => \$sample_id },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
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

sub add_program_metafile_to_sample_info {

## Function : Adds path and/or directory and/or file and/or version from programs to sample_info to track all metafiles and extract downstream
## Returns  :
## Arguments: $infile           => Infile for data at sample level {Optional}
##          : $directory        => Directory of the file
##          : $file             => File name
##          : $metafile_tag     => Id tag of meta file
##          : $path             => Path of file
##          : $processed_by     => Processed by
##          : $program_name     => Program name
##          : $version          => Version of file
##          : $sample_id        => Sample_id for data at sample level {Optional}
##          : $sample_info_href => Records on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $file;
    my $metafile_tag;
    my $infile;
    my $path;
    my $processed_by;
    my $program_name;
    my $sample_id;
    my $sample_info_href;
    my $version;

    my $tmpl = {
        directory => {
            strict_type => 1,
            store       => \$directory,
        },
        file => {
            strict_type => 1,
            store       => \$file,
        },
        infile       => { strict_type => 1, store => \$infile, },
        metafile_tag => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$metafile_tag,
        },
        path => {
            strict_type => 1,
            store       => \$path,
        },
        processed_by => {
            strict_type => 1,
            store       => \$processed_by,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_id        => { strict_type => 1, store => \$sample_id },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        version => {
            strict_type => 1,
            store       => \$version,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set the key and value pair to add to sample_info hash
    my %parameter = (
        directory    => $directory,
        file         => $file,
        path         => $path,
        processed_by => $processed_by,
        version      => $version,
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

1;
