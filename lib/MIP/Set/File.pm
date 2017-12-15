package MIP::Set::File;

use Carp;
use charnames qw{ :full :short };
use Cwd qw(abs_path);
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use Params::Check qw{ check allow last_error };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any };
use Readonly;

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ collect_path_entries set_absolute_path set_file_suffix set_merged_infile_prefix };
}

## Constants
Readonly my $NEWLINE => qq{\n};

sub collect_path_entries {

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

            collect_path_entries(
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

sub set_absolute_path {

## Function : Find aboslute path for supplied path or croaks and exists if path does not exists
## Returns  : $path (absolute path)
## Arguments: $parameter_name => Parameter to be evaluated
##          : $path           => Supplied path to be updated/evaluated

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $parameter_name;
    my $path;

    my $tmpl = {
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## For broadcasting later
    my $original_path = $path;

    ## Reformat to absolute path
    $path = abs_path($path);

    ## Something went wrong
    if ( not defined $path ) {

        croak(  q{Could not find absolute path for }
              . $parameter_name . q{: }
              . $original_path
              . q{. Please check the supplied path!} );
    }
    return $path;
}

sub set_file_suffix {

## Function : Set the current file suffix for this job id chain
## Returns  : $file_suffix
## Arguments: $file_suffix    => File suffix
##          : $job_id_chain   => Job id chain for program
##          : $parameter_href => Holds all parameters
##          : $suffix_key     => Suffix key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_suffix;
    my $job_id_chain;
    my $parameter_href;
    my $suffix_key;

    my $tmpl = {
        file_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$file_suffix,
            strict_type => 1,
        },
        job_id_chain => {
            defined     => 1,
            required    => 1,
            store       => \$job_id_chain,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        suffix_key => {
            defined     => 1,
            required    => 1,
            store       => \$suffix_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $parameter_href->{$suffix_key}{$job_id_chain} = $file_suffix;

    return $file_suffix;
}

sub set_merged_infile_prefix {

## Function : Set the merged infile prefix for sample id
## Returns  :
## Arguments: $file_info_href       => File info hash {REF}
##          : $merged_infile_prefix => Merged infile prefix
##          : $sample_id            => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $merged_infile_prefix;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        merged_infile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$merged_infile_prefix,
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

    $file_info_href->{$sample_id}{merged_infile} = $merged_infile_prefix;

    return;
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
