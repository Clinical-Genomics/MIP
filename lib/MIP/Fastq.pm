package MIP::Fastq;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use List::Util qw { any sum };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $CLOSE_PARENTHESIS $COMMA $EMPTY_STR $LOG_NAME $OPEN_PARENTHESIS $PIPE $PLUS $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      casava_header_features
      check_interleaved
      define_mip_fastq_file_features
      get_read_length
      get_stream_fastq_file_cmd
      parse_fastq_file_header_attributes
      parse_fastq_infiles
      parse_fastq_infiles_format };
}

## Constants
Readonly my $SUM_FOR_INTERLEAVED_DIRECTIONS => 3;

sub casava_header_features {

## Function : Define casava fastq file header formats features
## Returns  : %casava_header_feature
## Arguments:

    my ($arg_href) = @_;

    my %casava_header_feature = (
        q{1.4} => [qw{ instrument_id run_number lane tile x_pos y_pos index direction }],
        q{1.8} => [
            qw{ instrument_id run_number flowcell lane tile x_pos y_pos direction filtered control_bit index }
        ],
    );

    return %casava_header_feature;
}

sub check_interleaved {

## Function : Detect if fastq file is interleaved
## Returns  : "1(=interleaved)"
## Arguments: $file_path         => File to parse
##          : $read_file_command => Command used to read file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $read_file_command;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        read_file_command => {
            defined     => 1,
            required    => 1,
            store       => \$read_file_command,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Select relevant regexps to use
    my @regexps = qw{ get_fastq_header_v1.8_interleaved get_fastq_header_v1.4_interleaved };

    ## Store return from regexp
    my $fastq_read_direction;

  REGEXP:
    foreach my $regexp (@regexps) {

        ## Build regexp to find header features
        my @perl_commands = perl_nae_oneliners(
            {
                oneliner_name => $regexp,
            }
        );
        my $fastq_info_headers_cmd = qq{$read_file_command $file_path | @perl_commands;};

        my %return = child_process(
            {
                commands_ref => [$fastq_info_headers_cmd],
                process_type => q{ipc_cmd_run},
            }
        );

        $fastq_read_direction = $return{stdouts_ref}[0];
        last REGEXP if ($fastq_read_direction);
    }

    if ( not $fastq_read_direction ) {

        $log->fatal( q{Malformed fastq file: } . $file_path );
        $log->fatal(q{Could not find a read direction });
        exit 1;
    }

    my @fastq_read_directions = split //sxm, $fastq_read_direction;

    if ( any { /[^123]/sxm } @fastq_read_directions ) {

        $log->fatal(q{Malformed fastq file!});
        $log->fatal(
            q{Read direction is: } . join q{ and },
            @fastq_read_directions
              . q{, allowed entries are '1', '2', '3'. Please check fastq file}
              . $file_path
        );
        exit 1;
    }

    if ( sum(@fastq_read_directions) == $SUM_FOR_INTERLEAVED_DIRECTIONS ) {

        $log->info( q{Found interleaved fastq file: } . $file_path );
        return 1;
    }
    return;
}

sub define_mip_fastq_file_features {

## Function : Define MIP internal fastq file features derived from supplied fastq infile name
## Returns  : $mip_file_format, $mip_file_format_with_direction, $original_file_name_prefix, $run_barcode
## Arguments: $date               => Flow-cell sequencing date
##          : $direction          => Sequencing read direction
##          : $flowcell           => Flow-cell id
##          : $index              => The DNA library preparation molecular barcode
##          : $lane               => Flow-cell lane
##          : $original_file_name => Original file name
##          : $sample_id          => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $date;
    my $direction;
    my $flowcell;
    my $index;
    my $lane;
    my $original_file_name;
    my $sample_id;

    my $tmpl = {
        date => {
            defined     => 1,
            required    => 1,
            store       => \$date,
            strict_type => 1,
        },
        direction => {
            allow       => [ 1, 2 ],
            defined     => 1,
            required    => 1,
            store       => \$direction,
            strict_type => 1,
        },
        flowcell => {
            defined     => 1,
            required    => 1,
            store       => \$flowcell,
            strict_type => 1,
        },
        index => { defined => 1, required => 1, store => \$index, strict_type => 1, },
        lane  => {
            allow       => qr{ \A\d+\z }xsm,
            defined     => 1,
            required    => 1,
            store       => \$lane,
            strict_type => 1,
        },
        original_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$original_file_name,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $mip_file_format = join $UNDERSCORE,
      ( $sample_id, $date, $flowcell, $index, q{lane} . $lane );

    my $mip_file_format_with_direction = join $UNDERSCORE, ( $mip_file_format, $direction );

    my $original_file_name_prefix = substr $original_file_name, 0,
      index $original_file_name, q{.fastq};

    my $run_barcode = join $UNDERSCORE, ( $date, $flowcell, $lane, $index );

    return $mip_file_format, $mip_file_format_with_direction,
      $original_file_name_prefix, $run_barcode;
}

sub parse_fastq_file_header_attributes {

## Function : Get run info from fastq file header
## Returns  : @fastq_info_headers
## Arguments: $file_info_href    => File info hash {REF}
##          : $file_name         => Fast file name
##          : $file_path         => File path to parse
##          : $read_file_command => Command used to read file
##          : $sample_id         => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $file_name;
    my $file_path;
    my $read_file_command;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        read_file_command => {
            defined     => 1,
            required    => 1,
            store       => \$read_file_command,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Fastq qw{ casava_header_features };
    use MIP::File_info qw{ set_sample_file_attribute };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %casava_header_regexp = casava_header_features();

    my $fastq_info_header_string;

    my %fastq_header_info;

    my %regexp;

    ## Select relevant regexps to use
    @regexp{qw{ 1.4 1.8 }} = qw{ get_fastq_header_v1.4 get_fastq_header_v1.8 };

  REGEXP:
    while ( my ( $casava_version, $regexp ) = each %regexp ) {

        ## Build regexp to find header features
        my @perl_commands = perl_nae_oneliners(
            {
                oneliner_name => $regexp,
            }
        );

        ## Define cmd
        my $get_header_cmd = qq{$read_file_command $file_path | @perl_commands;};

        ## Collect fastq header info
        my %process_return = child_process(
            {
                commands_ref => [$get_header_cmd],
                process_type => q{ipc_cmd_run},
            }
        );

        $fastq_info_header_string = $process_return{stdouts_ref}[0];

        ## If successful regexp
        if ($fastq_info_header_string) {

            # Get header attributes
            my @attributes = @{ $casava_header_regexp{$casava_version} };

            # Split header string into array
            my @fastq_info_headers = split $SPACE, $fastq_info_header_string;

            # Add to hash to be returned
            @fastq_header_info{@attributes} = @fastq_info_headers;
            last REGEXP;
        }
    }

    if ( not $fastq_info_header_string ) {

        $log->fatal( q{Error parsing file header: } . $file_path );
        $log->fatal(q{Could not detect required sample sequencing run info from fastq file header});
        $log->fatal(q{Please proved MIP file in MIP file convention format to proceed});
        exit 1;
    }

    $log->warn(
q{Will add fake date '00010101' to follow file convention since this is not recorded in fastq header}
    );
    $fastq_header_info{date} = q{000101};

    if ( not exists $fastq_header_info{index} ) {

        # Special case since index is not present in fast headers casaava 1.4
        $fastq_header_info{index} = $EMPTY_STR;
    }

## Transfer to file_info hash
  ATTRIBUTE:
    while ( my ( $attribute, $attribute_value ) = each %fastq_header_info ) {

        set_sample_file_attribute(
            {
                attribute       => $attribute,
                attribute_value => $attribute_value,
                file_info_href  => $file_info_href,
                file_name       => $file_name,
                sample_id       => $sample_id,
            }
        );
    }
    ## Adds information derived from infile name to hashes
    $log->info( q{Found following information from fastq header: lane=}
          . $fastq_header_info{lane}
          . q{ flow-cell=}
          . $fastq_header_info{flowcell}
          . q{ index=}
          . $fastq_header_info{index}
          . q{ direction=}
          . $fastq_header_info{direction}, );
    return;
}

sub get_read_length {

## Function : Collect read length from a fastq infile
## Returns  : $read_length
## Arguments: $file_path => File to parse
##          : $read_file => Command used to read file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $read_file_command;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        read_file_command => {
            defined     => 1,
            required    => 1,
            store       => \$read_file_command,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Build regexp to find read length
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_name => q{get_fastq_read_length},
        }
    );

    my $read_length_cmd = qq{$read_file_command $file_path | @perl_commands};

    my %process_return = child_process(
        {
            commands_ref => [$read_length_cmd],
            process_type => q{ipc_cmd_run},
        }
    );

    ## Return read length
    return $process_return{stdouts_ref}[0];
}

sub get_stream_fastq_file_cmd {

## Function : Build command for streaming of chunk from fastq file(s)
## Returns  : @read_files_chunk_cmds
## Arguments: $fastq_files_ref => Fastq files {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fastq_files_ref;

    my $tmpl = {
        fastq_files_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$fastq_files_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Coreutils qw{ gnu_head gnu_tail };

    ## Constants
    Readonly my $BYTE_START_POS => 10_000;
    Readonly my $BYTE_STOP_POS  => $BYTE_START_POS * 300;

    my @read_files_chunk_cmds;

  FILE:
    foreach my $file_path ( @{$fastq_files_ref} ) {

        ## Check gzipped status of file path to choose correct cat binary (cat or gzip)
        ## Also prepend stream character
        my @cmds = _get_file_read_commands(
            {
                file_path => $file_path,
            }
        );

        push @cmds, $PIPE;

        ## Start of byte chunk
        push @cmds, gnu_tail( { number => $PLUS . $BYTE_START_POS, } );

        push @cmds, $PIPE;

        ## End of byte chunk
        push @cmds, gnu_head( { number => $BYTE_STOP_POS, } );

        $cmds[-1] = $cmds[-1] . $CLOSE_PARENTHESIS;

        ## Join for command line
        push @read_files_chunk_cmds, join $SPACE, @cmds;
    }
    return @read_files_chunk_cmds;
}

sub parse_fastq_infiles {

## Function : Parse fastq infiles for MIP processing
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href        => File info hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Fastq qw{ parse_fastq_file_header_attributes };
    use MIP::File_info qw{
      get_sample_file_attribute
      parse_files_compression_status
      parse_sample_fastq_file_attributes
    };
    use MIP::Sample_info qw{ set_infile_info };

    my $log = Log::Log4perl->get_logger($LOG_NAME);

  SAMPLE_ID:
    for my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        # Needed to be able to track when lanes are finished
        my $lane_tracker = 0;

        my %file_info_sample = get_sample_file_attribute(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );
        ## Unpack
        my $infiles_dir = $file_info_sample{mip_infiles_dir};

      INFILE:
        foreach my $file_name ( @{ $file_info_sample{mip_infiles} } ) {

            my %infile_info = parse_sample_fastq_file_attributes(
                {
                    file_info_href => $file_info_href,
                    file_name      => $file_name,
                    infiles_dir    => $infiles_dir,
                    sample_id      => $sample_id,
                }
            );

            ## If filename convention is followed
            if ( not exists $infile_info{date} ) {

                ## No regexp match i.e. file does not follow filename convention
                $log->warn(q{Will try to find mandatory information from fastq header});

                ## Get run info from fastq file header
                parse_fastq_file_header_attributes(
                    {
                        file_info_href    => $file_info_href,
                        file_name         => $file_name,
                        file_path         => catfile( $infiles_dir, $file_name ),
                        read_file_command => $infile_info{read_file_command},
                        sample_id         => $sample_id,
                    }
                );
            }

            ## Adds information derived from infile name to hashes
            $lane_tracker = set_infile_info(
                {
                    file_info_href   => $file_info_href,
                    file_name        => $file_name,
                    lane_tracker     => $lane_tracker,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            parse_files_compression_status(
                {
                    file_info_href => $file_info_href,
                    sample_id      => $sample_id,
                }
            );
        }
    }
    return;
}

sub parse_fastq_infiles_format {

## Function : Parse infile according to MIP filename convention
## Returns  : %infile_info or undef
## Arguments: $file_name => File name
##          : $sample_id => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;
    my $sample_id;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Validate::Case qw{ check_infile_contain_sample_id };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    # Store fastq file name features
    my %fastq_file_name;
    my @missing_feature;

    ## Define MIP fastq file name formats matching regexp
    my %fastq_file_name_regexp = _fastq_file_name_regexp();

    # Parse fastq file name
    my @file_features = $file_name =~ /$fastq_file_name_regexp{regexp}/sxm;

  FEATURE:
    while ( my ( $index, $feature ) = each @{ $fastq_file_name_regexp{features} } ) {

        ## Return undef if not all expected features found
        if ( not $file_features[$index] ) {
            push @missing_feature, $feature;
        }

        # Store feature
        $fastq_file_name{$feature} = $file_features[$index];
    }
    if (@missing_feature) {
        $log->warn(qq{Could not detect MIP file name convention for file: $file_name });
        $log->warn( q{Missing file name feature: } . join $COMMA . $SPACE, @missing_feature );

        ## Check that file name at least contains sample id
        return if ( $file_name =~ /$sample_id/sxm );

        $log->fatal(
            qq{Please check that the file name: $file_name contains the sample_id: $sample_id});
        exit 1;
    }
    ## Check that the sample_id provided and sample_id in infile name match
    check_infile_contain_sample_id(
        {
            infile_name      => $file_name,
            infile_sample_id => $fastq_file_name{infile_sample_id},
            sample_id        => $sample_id,
        }
    );
    return %fastq_file_name;
}

sub _fastq_file_name_regexp {

## Function : Define MIP fastq file name formats matching regexp
## Returns  : %fastq_file_name_regexp
## Arguments:

    my ($arg_href) = @_;

    my %fastq_file_name_regexp = (
        features => [qw{ lane date flowcell infile_sample_id index direction }],
        regexp   => q?(\d+)_(\d+)_([^_]+)_([^_]+)_([^_]+)_(\d).fastq?,
    );

    return %fastq_file_name_regexp;
}

sub _get_file_read_commands {

## Function : Check gzipped status of file path to choose correct cat binary (cat or gzip). Also prepend stream character.
## Returns  : @read_cmds
## Arguments: $file_path => Fastq File path to check status for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Coreutils qw{ gnu_cat };
    use MIP::Program::Gzip qw{ gzip };
    use MIP::Validate::Data qw{ %constraint };

    my @read_cmds;

    if ( $constraint{is_gzipped}->($file_path) ) {
        @read_cmds = gzip(
            {
                decompress       => 1,
                infile_paths_ref => [$file_path],
                stdout           => 1,
            }
        );
    }
    else {

        @read_cmds = gnu_cat( { infile_paths_ref => [$file_path], } );
    }

    $read_cmds[0] = q{<} . $OPEN_PARENTHESIS . $read_cmds[0];

    return @read_cmds;
}

1;
