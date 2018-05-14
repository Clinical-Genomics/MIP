package MIP::Get::File;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;
use List::MoreUtils qw{ any };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ get_exom_target_bed_file get_file_suffix get_merged_infile_prefix get_path_entries get_select_file_contigs get_seq_dict_contigs };
}

## Constants
Readonly my $COMMA   => q{,};
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub get_exom_target_bed_file {

## Function : Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash if supplied
## Returns  : $exome_target_bed_file
## Arguments: $exome_target_bed_href => Exome target bed files lnked to sample ids
##          : $file_ending           => File ending to add to file
##          : $log                   => Log object
##          : $sample_id             => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_href;
    my $file_ending;
    my $log;
    my $sample_id;

    my $tmpl = {
        exome_target_bed_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$exome_target_bed_href,
        },
        file_ending => { store => \$file_ending },
        log         => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## To track sample_ids with capture kits
    my %seen;

    ## Make local hash copy to keep in scope
    my %exome_target_bed = %{$exome_target_bed_href};

  BED_FILE:
    while ( my ( $exome_target_bed_file, $sample_id_string ) =
        each %exome_target_bed )
    {

        my @capture_kit_samples = split $COMMA, $sample_id_string;

        ## Count number of times sample_id has been seen
      SAMPLE:
        foreach my $samples (@capture_kit_samples) {

            $seen{$samples}++;
        }

        ## If capture_kit sample_id is associated with exome_target_bed file
        if ( any { $_ eq $sample_id } @capture_kit_samples ) {

            if ( defined $file_ending ) {

                $exome_target_bed_file .= $file_ending;
            }
            return $exome_target_bed_file;
        }
    }
    if ( not defined $seen{$sample_id} ) {

        $log->fatal(
            q{Could not detect }
              . $sample_id
              . q{ in '-exome_target_bed' associated files in sub routine get_exom_target_bed_file},
            $NEWLINE
        );
        exit 1;
    }
    return;
}

sub get_file_suffix {

## Function : Return the current file suffix for this jobid chain or program
## Returns  : $file_suffix
## Arguments: $job_id_chain   => Job id chain for program
##          : $parameter_href => Holds all parameters
##          : $program_name   => Program name
##          : $suffix_key     => Suffix key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $job_id_chain;
    my $parameter_href;
    my $program_name;
    my $suffix_key;

    my $tmpl = {
        jobid_chain    => { strict_type => 1, store => \$job_id_chain },
        program_name   => { strict_type => 1, store => \$program_name },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        suffix_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$suffix_key
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $file_suffix;

    ## Jobid chain specific suffix
    if ( defined $job_id_chain ) {

        $file_suffix = $parameter_href->{$suffix_key}{$job_id_chain};
    }
    elsif ( defined $program_name ) {
        ## Program  specific

        $file_suffix = $parameter_href->{$program_name}{$suffix_key};
    }

    ## If suffix was found
    if ( defined $file_suffix && $file_suffix ) {

        return $file_suffix;
    }
    else {
        ## Broadcast no suffix was found
        if ( defined $job_id_chain ) {

            say $log->fatal(
                q{Could not get requested infile_suffix for jobid_chain:}
                  . $job_id_chain );
        }
        elsif ( defined $program_name ) {

            say $log->fatal(
                q{Could not get requested infile_suffix for program:}
                  . $program_name );
        }
        exit 1;
    }
    return;
}

sub get_merged_infile_prefix {

## Function : Get the merged infile prefix for sample id
## Returns  : $merged_infile_prefix
## Arguments: $file_info_href => File info hash {REF}
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return $file_info_href->{$sample_id}{merged_infile};
}

sub get_path_entries {

## Function  : Collects all programs outfile path(s) created by MIP as Path->value located in %sample_info.
## Returns   :
## Arguments : $paths_ref        => Holds the collected paths {REF}
##           : $sample_info_href => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $paths_ref;
    my $sample_info_href;

    my $tmpl = {
        paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$paths_ref,
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

    ## Copy hash to enable recursive removal of keys
    my %info = %{$sample_info_href};

    ## Temporary array for collecting outDirectories within the same program
    my @outdirectories;

    ## Temporary array for collecting outfile within the same program
    my @outfiles;

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %info ) {

        if ( ref $value eq q{HASH} ) {

            get_path_entries(
                {
                    paths_ref        => $paths_ref,
                    sample_info_href => $value,
                }
            );
        }
        else {

            ## Required for first dry-run
            next KEY_VALUE_PAIR if ( not $value );

            ## Check if key is "path" and adds value to @paths_ref if true.
            _check_and_add_to_array(
                {
                    key       => $key,
                    paths_ref => $paths_ref,
                    value     => $value,
                }
            );

            ## Check if key is "outdirectory" or "outfile"  and adds joined value to @paths_ref if true.
            _collect_outfile(
                {
                    key                => $key,
                    paths_ref          => $paths_ref,
                    outdirectories_ref => \@outdirectories,
                    outfiles_ref       => \@outfiles,
                    value              => $value,
                }
            );

            delete $info{$value};
        }
    }
    return;
}

sub get_select_file_contigs {

## Function : Collects sequences contigs used in select file
## Returns  :
## Arguments: $log              => Log object
##          : $select_file_path => Select file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $select_file_path;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        select_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$select_file_path
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    my $find_contig_name = q?perl -nae ?;

    ## Get contig name
    $find_contig_name .= q?'if ($_=~/ contig=(\w+) /xsm) { ?;

    ## Alias capture
    $find_contig_name .= q?my $contig_name = $1; ?;

    ## Write contig name and comma
    $find_contig_name .= q?print $contig_name, q{,};} ?;

    ## Quit if #CHROM found in line
    $find_contig_name .= q?if($_=~/ [#]CHROM /xsm) {last;}' ?;

    ## Returns a comma seperated string of sequence contigs from file
    my @contigs = `$find_contig_name $select_file_path `;

    @contigs = split $COMMA, join $COMMA, @contigs;

    if ( not @contigs ) {

        $log->fatal(
q{Could not detect any '##contig' in meta data header in select file: }
              . $select_file_path );
        exit 1;
    }
    return @contigs;
}

sub get_seq_dict_contigs {

## Function : Collects sequences contigs used in analysis from human genome sequence dictionnary associated with $human_genome_reference
## Returns  :
## Arguments: $dict_file_path => Dict file path
##          : $log            => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dict_file_path;
    my $log;

    my $tmpl = {
        dict_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$dict_file_path,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### Build regexp to find contig names
    ## Execute perl
    my $find_contig_name = q?perl -nae '?;

    ## Find contig line
    $find_contig_name .= q?if($F[0]=~/^\@SQ/) { ?;

    ## Collect contig name
    $find_contig_name .= q? if($F[1]=~/SN\:(\S+)/) { ?;

    ## Alias capture
    $find_contig_name .= q?my $contig_name = $1; ?;

    ## Write to STDOUT
    $find_contig_name .= q?print $contig_name, q{,};} }' ?;

    ## Returns a comma seperated string of sequence contigs from dict file
    my @contigs = `$find_contig_name $dict_file_path `;

    ## Save contigs
    @contigs = split $COMMA, join $COMMA, @contigs;

    if ( not @contigs ) {

        log->fatal( q{Could not detect any 'SN:contig_names' in dict file: }
              . $dict_file_path );
        exit 1;
    }
    return @contigs;
}

sub _check_and_add_to_array {

## Function  : Check if Key name is "path" and adds to @paths_ref if true.
## Returns   :
## Arguments : $keyName   => Hash key
##           : $paths_ref => Holds the collected paths {REF}
##           : $value     => Hash value

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $key;
    my $paths_ref;
    my $value;

    my $tmpl = {
        key =>
          { defined => 1, required => 1, store => \$key, strict_type => 1, },
        paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$paths_ref,
            strict_type => 1,
        },
        value => { required => 1, store => \$value, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( $key ne q{path} );

    ## Do not add same path twice
    if ( not any { $_ eq $value } @{$paths_ref} ) {

        push @{$paths_ref}, $value;
    }
    return;
}

sub _collect_outfile {

## Function  : Check if Key name is "outdirectory" or "outfile"  and adds to @paths_ref if true.
## Returns   :
## Arguments : $key                => Hash key
##           : $outdirectories_ref => Holds temporary outdirectory path(s) {Optional, REF}
##           : $outfiles_ref       => Holds temporary outdirectory path(s) {Optional, REF}
##           : $value              => Hash value
##           : $paths_ref          => Holds the collected paths {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $key;
    my $outdirectories_ref;
    my $outfiles_ref;
    my $paths_ref;
    my $value;

    my $tmpl = {
        paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$paths_ref,
            strict_type => 1,
        },
        outdirectories_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$outdirectories_ref,
            strict_type => 1,
        },
        outfiles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$outfiles_ref,
            strict_type => 1,
        },
        value => { defined => 1, required => 1, store => \$value, },
        key =>
          { defined => 1, store => \$key, required => 1, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $key eq q{outdirectory} ) {

        push @{$outdirectories_ref}, $value;
    }
    if ( $key eq q{outfile} ) {

        push @{$outfiles_ref}, $value;
    }

    ## Both outdirectory and outfile have been collected, time to join
    if ( @{$outdirectories_ref} && @{$outfiles_ref} ) {

        my $path = catfile( $outdirectories_ref->[0], $outfiles_ref->[0] );

        ## Do not add same path twice
        if ( not any { $_ eq $path } @{$paths_ref} ) {

            push @{$paths_ref},
              catfile( $outdirectories_ref->[0], $outfiles_ref->[0] );

            ## Restart
            @{$outdirectories_ref} = ();
            @{$outfiles_ref}       = ();
        }
    }
    return;
}

1;
