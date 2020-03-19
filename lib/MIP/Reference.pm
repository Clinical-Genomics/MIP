package MIP::Reference;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catfile };
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
use MIP::Constants qw{ $COMMA $LOG_NAME $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_human_genome_file_endings
      get_dict_contigs
      parse_meta_file_suffixes
      update_exome_target_bed
      write_contigs_size_file
    };
}

sub check_human_genome_file_endings {

## Function : Check the existance of associated human genome files
## Returns  :
## Arguments: $human_genome_reference_file_endings_ref => Human genome reference file endings
##          : $human_genome_reference_path             => Human genome reference file path
##          : $parameter_href                          => Parameter hash {REF}
##          : $parameter_name                          => The parameter under evaluation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $human_genome_reference_file_endings_ref;
    my $human_genome_reference_path;
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        human_genome_reference_file_endings_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_file_endings_ref,
            strict_type => 1,
        },
        human_genome_reference_path => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_path,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => { store => \$parameter_name, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_existance };
    use MIP::Parameter qw{ set_parameter_build_file_status };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Count the number of files that exists
    my $existence_check_counter = 0;

  FILE_ENDING:
    foreach my $file_ending ( @{$human_genome_reference_file_endings_ref} ) {

        ## Duplicate to keep fresh with each loop
        my $path = $human_genome_reference_path;

        ## Dict requires no fastq(.gz) ending
        if ( $file_ending eq q{.dict} ) {

            ## Removes ".file_ending" in filename.FILENDING(.gz)
            my ( $file, $dir_path ) =
              fileparse( $path, qr/ [.]fasta | [.]fasta[.]gz /sxm );
            $path = catfile( $dir_path, $file );
        }

        ## Add current ending
        $path = $path . $file_ending;

        my ($does_exist) = check_filesystem_objects_existance(
            {
                object_name    => $path,
                object_type    => q{file},
                parameter_name => $parameter_name,
            }
        );

        ## Sum up the number of file that exists
        $existence_check_counter = $existence_check_counter + $does_exist;
    }
    ## Files need to be built
    if ( $existence_check_counter != scalar @{$human_genome_reference_file_endings_ref} )
    {

        set_parameter_build_file_status {
            (
                parameter_href => $parameter_href,
                parameter_name => $parameter_name,
                status         => 1,
            )
        };
        return;
    }

    # All files exist in this check
    set_parameter_build_file_status {
        (
            parameter_href => $parameter_href,
            parameter_name => $parameter_name,
            status         => 0,
        )
    };
    return;
}

sub get_dict_contigs {

## Function : Collects sequence contigs used in analysis from human genome sequence
##          : dictionnary (.dict file)
## Returns  : @contigs
## Arguments: $dict_file_path => Dict file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dict_file_path;

    my $tmpl = {
        dict_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$dict_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Build regexp to find contig names
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_name  => q{get_dict_contigs},
            stdinfile_path => $dict_file_path,
        }
    );

    my @get_dict_contigs_cmds = join $SPACE, ( @perl_commands, );

    # System call
    my %return = child_process(
        {
            commands_ref => \@get_dict_contigs_cmds,
            process_type => q{ipc_cmd_run},
        }
    );

    # Save contigs
    my @contigs = split $COMMA, join $COMMA, @{ $return{stdouts_ref} };

    #my @contigs = split $COMMA, join $COMMA, @{$stdout_buf_ref};

    return @contigs if (@contigs);

    $log->fatal(
        q{Could not detect any 'SN:contig_names' in dict file: } . $dict_file_path );
    exit 1;
}

sub parse_meta_file_suffixes {

## Function : Checks files to be built by combining object name prefix with suffix.
## Returns  :
## Arguments: $active_parameter_href  => Active parameters for this analysis hash {REF}
##          : $file_name              => File name
##          : $meta_file_suffixes_ref => Reference to the meta file suffixes to be added to the file name {REF}
##          : $parameter_href         => Parameter hash {REF}
##          : $parameter_name         => MIP parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_name;
    my $meta_file_suffixes_ref;
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        meta_file_suffixes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$meta_file_suffixes_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_existance };
    use MIP::Parameter qw{ get_parameter_attribute set_parameter_build_file_status };

    ## Count the number of files that exists
    my $existence_check_counter = 0;

    my $build_status = 0;

    ## Get parameter object type i.e file or directory
    my $object_type = get_parameter_attribute(
        {
            attribute      => q{exists_check},
            parameter_href => $parameter_href,
            parameter_name => $parameter_name,
        }
    );

  FILE_SUFFIX:
    foreach my $file_suffix ( @{$meta_file_suffixes_ref} ) {

        my ($exist) = check_filesystem_objects_existance(
            {
                object_name    => catfile( $file_name . $file_suffix ),
                object_type    => $object_type,
                parameter_name => $parameter_name,
            }
        );
        ## Sum up the number of file that exists
        $existence_check_counter = $existence_check_counter + $exist;
    }

    ## Files need to be built
    if ( $existence_check_counter != scalar @{$meta_file_suffixes_ref} ) {

        $build_status = 1;
    }

    # Set build status for parameter
    set_parameter_build_file_status(
        {
            parameter_href => $parameter_href,
            parameter_name => $parameter_name,
            status         => $build_status,
        }
    );
    return;
}

sub update_exome_target_bed {

## Function : Update exome_target_bed files with human genome reference source and version
## Returns  :
## Arguments: $exome_target_bed_file_href     => Exome target bed {REF}
##          : $human_genome_reference_source  => Human genome reference source
##          : $human_genome_reference_version => Human genome reference version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_file_href;
    my $human_genome_reference_source;
    my $human_genome_reference_version;

    my $tmpl = {
        exome_target_bed_file_href =>
          { required => 1, store => \$exome_target_bed_file_href, },
        human_genome_reference_source => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_source,
            strict_type => 1,
        },
        human_genome_reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  EXOME_FILE:
    foreach my $exome_target_bed_file ( keys %{$exome_target_bed_file_href} ) {

        my $original_file_name = $exome_target_bed_file;

        ## Replace with actual version
        if ( $exome_target_bed_file =~
            s/genome_reference_source/$human_genome_reference_source/xsm
            && $exome_target_bed_file =~ s/_version/$human_genome_reference_version/xsm )
        {

            ## The delete operator returns the value being deleted
            ## i.e. updating hash key while preserving original info
            $exome_target_bed_file_href->{$exome_target_bed_file} =
              delete $exome_target_bed_file_href->{$original_file_name};
        }
    }
    return;
}

sub write_contigs_size_file {

## Function : Write contig size file from human genome sequence fai (.fai) file
## Returns  :
## Arguments: $fai_file_path => Fai file path
##          : $outfile_path  => Chromosome size file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fai_file_path;
    my $outfile_path;

    my $tmpl = {
        fai_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$fai_file_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Build regexp to find contig names
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_name   => q{write_contigs_size_file},
            stdinfile_path  => $fai_file_path,
            stdoutfile_path => $outfile_path,
        }
    );

    my @write_contigs_size_cmd = join $SPACE, ( @perl_commands, );

    # System call
    my %process_return = child_process(
        {
            commands_ref => \@write_contigs_size_cmd,
            process_type => q{ipc_cmd_run},
            verbose      => 0,
        }
    );

    return if ( not @{ $process_return{stderrs_ref} } );

    $log->fatal(q{Could not write contigs size file});
    $log->fatal( q{Error: } . join $NEWLINE, @{ $process_return{stderrs_ref} } );
    exit 1;
}

1;
