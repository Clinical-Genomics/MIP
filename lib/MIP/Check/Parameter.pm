package MIP::Check::Parameter;

use Carp;
use charnames qw{ :full :short };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Email::Valid;
use Readonly;
use List::MoreUtils qw { any };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_allowed_array_values
      check_allowed_temp_directory
      check_aligner
      check_cmd_config_vs_definition_file
      check_email_address
      check_gzipped
      check_infile_contain_sample_id
      check_infiles
      check_parameter_hash
      check_program_exists_in_hash
      check_prioritize_variant_callers
      check_program_mode
      check_sample_ids
      check_sample_id_in_hash_parameter
      check_sample_id_in_hash_parameter_path
      check_snpsift_keys
      check_vep_directories
    };
}

## Constants
Readonly my $COMMA         => q{,};
Readonly my $DOLLAR_SIGN   => q{$};
Readonly my $DOT           => q{.};
Readonly my $FORWARD_SLASH => q{/};
Readonly my $NEWLINE       => qq{\n};
Readonly my $PIPE          => q{|};
Readonly my $SINGLE_QUOTE  => q{'};
Readonly my $SPACE         => q{ };
Readonly my $UNDERSCORE    => q{_};

sub check_allowed_array_values {

## Function : Check that the array values are allowed
## Returns  :
## Arguments: $allowed_values_ref => Allowed values for parameter {REF}
##          : $values_ref         => Values for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $allowed_values_ref;
    my $values_ref;

    my $tmpl = {
        allowed_values_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$allowed_values_ref,
            strict_type => 1,
        },
        values_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$values_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %is_allowed;

    # Remake allowed values into keys in is_allowed hash
    map { $is_allowed{$_} = undef } @{$allowed_values_ref};

  VALUES:
    foreach my $value ( @{$values_ref} ) {

        # Test if value is allowed
        if ( not exists $is_allowed{$value} ) {

            return 0;
        }
    }

    # All ok
    return 1;
}

sub check_allowed_temp_directory {

## Function : Check that the temp directory value is allowed
## Returns  :
## Arguments: $log            => Log object
##          : $temp_directory => Temp directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $temp_directory;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %is_not_allowed = (
        $FORWARD_SLASH . q{scratch}                  => undef,
        $FORWARD_SLASH . q{scratch} . $FORWARD_SLASH => undef,
    );

    # Test if value is allowed
    if ( exists $is_not_allowed{$temp_directory} ) {

        $log->fatal( qq{$SINGLE_QUOTE--temp_directory }
              . $temp_directory
              . qq{$SINGLE_QUOTE is not allowed because MIP will remove the temp directory after processing.}
        );
        exit 1;
    }

    # All ok
    return 1;
}

sub check_aligner {

## Function : Check that the correct number of aligners is used in MIP and sets the outaligner_dir flag accordingly.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##          : $log                   => Log object
##          : $outaligner_dir        => Outaligner_dir used in the analysis
##          : $parameter_href        => Parameter hash {REF}
##          : $verbose               => Verbosity level

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $log;
    my $parameter_href;

    ## Default(s)
    my $outaligner_dir;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        broadcasts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$broadcasts_ref,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        verbose => {
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %aligner;

  ALIGNER:
    foreach my $aligner ( @{ $parameter_href->{dynamic_parameter}{aligners} } )
    {

        ## Active aligner
        if ( $active_parameter_href->{$aligner} ) {

            # Increment aligner count
            $aligner{total_active_aligner_count}++;

            # Store aligner
            push @{ $aligner{active_aligners} }, $aligner;

            # Set active aligner for downstream use
            $parameter_href->{active_aligner} = $aligner;

            if ( not defined $outaligner_dir ) {

                # Set outaligner_dir parameter depending on active aligner
                $active_parameter_href->{outaligner_dir} = $outaligner_dir =
                  $parameter_href->{$aligner}{outdir_name};

                next ALIGNER if ( not $verbose );

                my $info = q{Set outaligner_dir to: } . $outaligner_dir;

                ## Add info to broadcasts
                push @{$broadcasts_ref}, $info;
            }
        }
    }

    if ( exists $aligner{total_active_aligner_count}
        and $aligner{total_active_aligner_count} > 1 )
    {

        $log->fatal( q{You have activate more than 1 aligner: }
              . join( q{, }, @{ $aligner{active_aligners} } )
              . q{. MIP currently only supports 1 aligner per analysis.} );
        exit 1;
    }
    return;
}

sub check_cmd_config_vs_definition_file {

## Function : Compare keys from config and cmd with definitions file
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @allowed_unique_keys =
      ( q{vcfparser_outfile_count}, $active_parameter_href->{family_id} );
    my @unique;

  ACTIVE_PARAMETER:
    foreach my $key ( keys %{$active_parameter_href} ) {

        ## Parameters from definitions file
        if ( not exists $parameter_href->{$key} ) {

            push @unique, $key;
        }
    }
  UNIQUE_KEYS:
    foreach my $unique_key (@unique) {

        ## Do not print if allowed_unique_keys that have been created dynamically from previous runs
        if ( not any { $_ eq $unique_key } @allowed_unique_keys ) {

            say {*STDERR} q{Found illegal key: }
              . $unique_key
              . q{ in config file or command line that is not defined in define_parameters.yaml};
            croak();
        }
    }
    return;
}

sub check_email_address {

## Function : Check the syntax of the email adress is valid and has a mail host.
## Returns  :
## Arguments: $email => The email adress
##          : $log   => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $email;
    my $log;

    my $tmpl = {
        email => {
            required    => 1,
            store       => \$email,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $email );

    ## Check syntax and mail host
    my $address = Email::Valid->address(
        -address => $email,
        -mxcheck => 1,
    );
    if ( not defined $address ) {

        $log->fatal( q{The supplied email: }
              . $email
              . q{ seem to be malformed according to }
              . Email::Valid->details() );
        exit 1;
    }
    return 1;
}

sub check_gzipped {

## Function : Check if a file is gzipped.
## Returns  : "0 (=uncompressed)| 1 (=compressed)"
## Arguments: $file_name => File name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $file_compression_status = 0;

    if ( $file_name =~ / [.]gz$ /xms ) {

        $file_compression_status = 1;
    }
    return $file_compression_status;
}

sub check_infile_contain_sample_id {

## Function : Check that the sample_id provided and sample_id in infile name match.
## Returns  :
## Arguments: $infile_name      => Infile name
##          : $infile_sample_id => Sample_id collect with regexp from infile
##          : $log              => Log object
##          : $sample_id        => Sample id from user
##          : $sample_ids_ref   => Sample ids from user

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_name;
    my $infile_sample_id;
    my $log;
    my $sample_id;
    my $sample_ids_ref;

    my $tmpl = {
        infile_name => {
            defined     => 1,
            required    => 1,
            store       => \$infile_name,
            strict_type => 1,
        },
        infile_sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$infile_sample_id,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Track seen sample ids
    my %seen;

    ## Increment for all sample ids
    map { $seen{$_}++ } ( @{$sample_ids_ref}, $infile_sample_id );

    if ( not $seen{$infile_sample_id} > 1 ) {

        $log->fatal( $sample_id
              . q{ supplied and sample_id }
              . $infile_sample_id
              . q{ found in file : }
              . $infile_name
              . q{ does not match. Please rename file to match sample_id: }
              . $sample_id );
        exit 1;
    }
    return 1;
}

sub check_infiles {

## Function : Check infiles found and that they contain sample_id
## Returns  :
## Arguments: $infiles_ref      => Infiles to check {REF}
##          : $infile_directory => Infile directory
##          : $sample_id        => Sample id
##          : $log              => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infiles_ref;
    my $infile_directory;
    my $log;
    my $sample_id;

    my $tmpl = {
        infile_directory => {
            defined     => 1,
            required    => 1,
            store       => \$infile_directory,
            strict_type => 1,
        },
        infiles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infiles_ref,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## No "*.fastq*" infiles
    if ( not @{$infiles_ref} ) {

        $log->fatal(
            q{Could not find any '.fastq' files in supplied infiles directory }
              . $infile_directory,
        );
        exit 1;
    }

    ## Check that infiledirs/infile contains sample_id in filename
  INFILE:
    foreach my $infile ( @{$infiles_ref} ) {

        next INFILE if ( $infile =~ /$sample_id/sxm );

        $log->fatal( q{Could not detect sample_id: }
              . $sample_id
              . q{ in supplied infile: }
              . catfile( $infile_directory, $infile ) );
        $log->fatal(
q{Check that: '--sample_ids' and '--infile_dirs' contain the same sample_id and that the filename of the infile contains the sample_id.},
        );
        exit 1;
    }
    return 1;
}

sub check_parameter_hash {

## Function : Evaluate parameters in parameters hash
## Returns  :
## Arguments: $file_path              => Path to yaml file
##          : $mandatory_key_href     => Hash with mandatory key {REF}
##          : $non_mandatory_key_href => Hash with non mandatory key {REF}
##          : $parameter_href         => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $mandatory_key_href;
    my $non_mandatory_key_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        mandatory_key_href => {
            default     => {},
            required    => 1,
            store       => \$mandatory_key_href,
            strict_type => 1,
        },
        non_mandatory_key_href => {
            default     => {},
            required    => 1,
            store       => \$non_mandatory_key_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check that mandatory keys exists for each parameter
    _check_parameter_mandatory_keys_exits(
        {
            file_path          => $file_path,
            mandatory_key_href => $mandatory_key_href,
            parameter_href     => $parameter_href,
        }
    );

    ## Test both mandatory and non_mandatory keys data type and values
    my @arguments = ( \%{$mandatory_key_href}, \%{$non_mandatory_key_href} );

  ARGUMENT_HASH_REF:
    foreach my $argument_href (@arguments) {

        ## Mandatory keys
        _check_parameter_keys(
            {
                file_path      => $file_path,
                key_href       => $argument_href,
                parameter_href => $parameter_href,
            }
        );
    }
    return 1;
}

sub check_program_exists_in_hash {

## Function : Test if parameter "program name" from query parameter exists truth hash
## Returns  :
## Arguments: $log            => Log object
##          : $parameter_name => Parameter name
##          : $query_ref      => Query (ARRAY|HASH|SCALAR) {REF}
##          : $truth_href     => Truth hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $parameter_name;
    my $query_ref;
    my $truth_href;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_name =>
          { defined => 1, required => 1, store => \$parameter_name, },
        truth_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$truth_href,
            strict_type => 1,
        },
        query_ref => {
            defined  => 1,
            required => 1,
            store    => \$query_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $error_msg =
      qq{$SINGLE_QUOTE - Does not exist as module program parameter in MIP};

    if ( ref $query_ref eq q{HASH} ) {

      PROGRAM_NAME:
        foreach my $program_name ( keys %{$query_ref} ) {

            next PROGRAM_NAME if ( exists $truth_href->{$program_name} );

            $log->fatal( $parameter_name
                  . qq{ key $SINGLE_QUOTE}
                  . $program_name
                  . $error_msg );
            exit 1;
        }
    }
    if ( ref $query_ref eq q{ARRAY} ) {

      PROGRAM_NAME:
        foreach my $program_name ( @{$query_ref} ) {

            next PROGRAM_NAME if ( exists $truth_href->{$program_name} );

            $log->fatal( $parameter_name
                  . qq{ element $SINGLE_QUOTE}
                  . $program_name
                  . $error_msg );
            exit 1;
        }
    }
    if ( ref $query_ref eq q{SCALAR} ) {

        return if ( exists $truth_href->{$parameter_name} );

        $log->fatal( $parameter_name
              . qq{ element $SINGLE_QUOTE}
              . $parameter_name
              . $error_msg );
        exit 1;
    }
    return;
}

sub check_prioritize_variant_callers {

## Function : Check that all active variant callers have a prioritization order and that the prioritization elements match a supported variant caller.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $parameter_href        => Parameter hash {REF}
##          : $parameter_name        => Parameter name
##          : $variant_callers_ref   => Variant callers to check {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $parameter_href;
    my $parameter_name;
    my $variant_callers_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
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
        variant_callers_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$variant_callers_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @priority_order_names =
      split $COMMA, $active_parameter_href->{$parameter_name};

    ## Alias
    my @variant_caller_aliases;

    ## Check that all active variant callers have a priority order
  CALLER:
    foreach my $variant_caller ( @{$variant_callers_ref} ) {

        my $variant_caller_alias =
          $parameter_href->{$variant_caller}{outdir_name};
        push @variant_caller_aliases, $variant_caller_alias;

        ## Only active programs
        if ( $active_parameter_href->{$variant_caller} ) {

            ## If variant caller alias is not part of priority order names
            if ( not any { $_ eq $variant_caller_alias } @priority_order_names )
            {

                $log->fatal( $parameter_name
                      . q{ does not contain active variant caller: '}
                      . $variant_caller_alias
                      . $SINGLE_QUOTE );
                exit 1;
            }
        }
        else {
            ## Only NOT active programs

            ## If variant caller alias is part of priority order names
            if ( any { $_ eq $variant_caller_alias } @priority_order_names ) {

                $log->fatal( $parameter_name
                      . q{ contains deactivated variant caller: '}
                      . $variant_caller_alias
                      . $SINGLE_QUOTE );
                exit 1;
            }
        }
    }

    ## Check that prioritize string contains valid variant call names
  PRIO_CALL:
    foreach my $prioritize_call (@priority_order_names) {

        # If priority order names is not part of variant caller alias
        if ( not any { $_ eq $prioritize_call } @variant_caller_aliases ) {

            $log->fatal( $parameter_name . q{: '}
                  . $prioritize_call
                  . q{' does not match any supported variant caller: '}
                  . join( $COMMA, @variant_caller_aliases )
                  . $SINGLE_QUOTE );
            exit 1;
        }
    }
    return 1;
}

sub check_program_mode {

## Function : Check correct value for program mode in MIP.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set allowed values
    my %is_allowed = map { $_ => 1 } ( 0 .. 2 );

  PROGRAM:
    foreach my $program ( @{ $parameter_href->{dynamic_parameter}{program} } ) {

        ## Alias
        my $program_mode = $active_parameter_href->{$program};

        next PROGRAM if ( $is_allowed{$program_mode} );

        #If not an allowed value in active parameters
        $log->fatal(
            $SINGLE_QUOTE
              . $active_parameter_href->{$program}
              . q{' Is not an allowed mode for program '--}
              . $program
              . q{'. Set to: }
              . join $PIPE,
            ( sort keys %is_allowed )
        );
        exit 1;
    }
    return 1;
}

sub check_sample_ids {

## Function : Test that the family_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
## Returns  :
## Arguments: $family_id      => Family id
##          : $log            => Log object
##          : $sample_ids_ref => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $family_id;
    my $log;
    my $sample_ids_ref;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        family_id => {
            defined     => 1,
            required    => 1,
            store       => \$family_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Hash to test duplicate sample_ids later
    my %seen;

    if ( not @{$sample_ids_ref} ) {

        $log->fatal(q{Please provide sample_id(s)});
        exit 1;
    }

  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        ## Increment instance to check duplicates later
        $seen{$sample_id}++;

        ## Family_id cannot be the same as sample_id
        if ( $family_id eq $sample_id ) {

            $log->fatal( q{Family_id: }
                  . $family_id
                  . q{ equals sample_id: }
                  . $sample_id
                  . q{. Please make sure that the family_id and sample_id(s) are unique.}
            );
            exit 1;
        }
        ## Check for unique sample_ids
        if ( $seen{$sample_id} > 1 ) {

            $log->fatal( q{Sample_id: } . $sample_id . q{ is not uniqe.} );
            exit 1;
        }
        ## Sample_id contains "_", not allowed in filename convention
        if ( $sample_id =~ /$UNDERSCORE/sxm ) {

            $log->fatal( q{Sample_id: }
                  . $sample_id
                  . q{ contains '_'. Please rename sample_id according to MIP's filename convention, removing the '_'.}
            );
            exit 1;
        }
    }
    return 1;
}

sub check_sample_id_in_hash_parameter {

## Function : Check sample_id provided in hash parameter is included in the analysis and only represented once
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $parameter_href        => Holds all parameters {REF}
##          : $parameter_names_ref   => Parameter name list {REF}
##          : $sample_ids_ref        => Array to loop in for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $parameter_names_ref;
    my $parameter_href;
    my $sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_names_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parameter_names_ref,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter_name ( @{$parameter_names_ref} ) {

        ## Skip undef parameters in current analysis
        next PARAMETER
          if ( not defined $active_parameter_href->{$parameter_name} );

      SAMPLE_ID:
        foreach my $sample_id ( @{$sample_ids_ref} ) {

            ## Unpack
            my $sample_id_value =
              $active_parameter_href->{$parameter_name}{$sample_id};
            my $is_mandatory = $parameter_href->{$parameter_name}{mandatory};

            ## Check that a value exists
            if ( not defined $sample_id_value ) {

                next PARAMETER
                  if ( defined $is_mandatory and $is_mandatory eq q{no} );

                $log->fatal(
                    q{Could not find value for }
                      . $sample_id
                      . q{ for parameter '--}
                      . $parameter_name
                      . $SINGLE_QUOTE
                      . q{. Provided sample_ids for parameter are: }
                      . join $COMMA
                      . $SPACE,
                    ( keys %{ $active_parameter_href->{$parameter_name} } )
                );
                exit 1;
            }
        }
    }
    return 1;
}

sub check_sample_id_in_hash_parameter_path {

## Function : Check sample_id provided in hash path parameter is included in the analysis and only represented once
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $parameter_names_ref   => Parameter name list {REF}
##          : $sample_ids_ref        => Array to loop in for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $parameter_names_ref;
    my $sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_names_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parameter_names_ref,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter_name ( @{$parameter_names_ref} ) {

        # Hash to test duplicate sample_ids later
        my %seen;

      SAMPLE_STR:
        foreach my $sample_id_str (
            keys %{ $active_parameter_href->{$parameter_name} } )
        {

            ## Get sample ids for parameter
            my @parameter_samples =
              split $COMMA,
              $active_parameter_href->{$parameter_name}{$sample_id_str};

          SAMPLE_ID:
            foreach my $sample_id (@parameter_samples) {

                # Increment instance to check duplicates later
                $seen{$sample_id}++;

                ## Check sample_id are unique
                if ( $seen{$sample_id} > 1 ) {

                    $log->fatal(
                        q{Sample_id: }
                          . $sample_id
                          . q{ is not uniqe in '--}
                          . $parameter_name . q{': }
                          . $sample_id_str . q{=}
                          . join $COMMA,
                        @parameter_samples,
                    );
                    exit 1;
                }
            }
        }
        ## Check all sample ids are present in parameter string
      SAMPLE_ID:
        foreach my $sample_id ( @{$sample_ids_ref} ) {

            ## If sample_id is not present in parameter_name hash
            if ( not any { $_ eq $sample_id } ( keys %seen ) ) {

                $log->fatal(
                    q{Could not detect }
                      . $sample_id
                      . q{ for '--}
                      . $parameter_name
                      . q{'. Provided sample_ids are: }
                      . join $COMMA
                      . $SPACE,
                    ( keys %seen ),
                );
                exit 1;
            }
        }
    }
    return 1;
}

sub check_snpsift_keys {

## Function : Check that the supplied snpsift outinfo keys match annotation files
## Returns  :
## Arguments: $log                                 => Log object
##          : $snpsift_annotation_files_href       => Snpsift annotation files {REF}
##          : $snpsift_annotation_outinfo_key_href => File and outinfo key to add to vcf {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $snpsift_annotation_files_href;
    my $snpsift_annotation_outinfo_key_href;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        snpsift_annotation_files_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$snpsift_annotation_files_href,
            strict_type => 1,
        },
        snpsift_annotation_outinfo_key_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$snpsift_annotation_outinfo_key_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  FILE:
    foreach my $file ( keys %{$snpsift_annotation_outinfo_key_href} ) {

        ## Matching files
        next FILE if ( exists $snpsift_annotation_files_href->{$file} );

        ## Else croak and exist
        $log->fatal( q{The supplied snpsift_annotation_outinfo_key file: }
              . $file
              . q{ does not match any file in '--snpsift_annotation_files'} );
        $log->fatal(
            q{Supplied snpsift_annotation_files files:}
              . $NEWLINE
              . join $NEWLINE,
            keys %{$snpsift_annotation_files_href}
        );
        exit 1;
    }
    return 1;
}

sub check_vep_directories {

## Function : Compare VEP directory and VEP chache versions. Exit if non-match
## Returns  :
## Arguments: $log                 => Log object
##          : $vep_directory_path  => VEP directory path {REF}
##          : $vep_directory_cache => VEP cache directory path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $vep_directory_path;
    my $vep_directory_cache;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        vep_directory_path => {
            defined     => 1,
            required    => 1,
            store       => \$vep_directory_path,
            strict_type => 1,
        },
        vep_directory_cache => {
            defined     => 1,
            required    => 1,
            store       => \$vep_directory_cache,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use File::Find::Rule;
    use Cwd qw{ abs_path };

    ## Get VEP version from file
    my $vep_version_file =
      catfile( $vep_directory_path, $DOT . q{version}, q{ensembl} );
    open my $vep_version_fh, q{<}, $vep_version_file;
    my $vep_version_line = <$vep_version_fh>;
    close $vep_version_fh;

    ## Capture version number
    my ($vep_version) = $vep_version_line =~ /( \d+ )/xms;

    ## Check that a version number was picked up
    if ( not $vep_version ) {
        $log->warn(
q{Could not retrieve VEP version. Skipping checking that VEP api and cache matches.}
        );
        return;
    }

    ## Get absolute path to VEP cache as it is commonly linked
    my $vep_cache_dir_path =
      abs_path( catdir( $vep_directory_cache, q{homo_sapiens} ) );

    ## Get folders in cache directory
    # Build rule
    my $rule = File::Find::Rule->new();

    # Find directories
    $rule->directory;

    # Only get directories in the current folder
    $rule->maxdepth(1);

    # Get relative paths
    $rule->relative;

    # Apply rule to find directories
    my @vep_cache_version_folders = $rule->in($vep_cache_dir_path);

    ## Check that
    if ( not @vep_cache_version_folders ) {
        $log->warn(
q{Could not retrieve VEP cache version. Skipping checking that VEP api and cache matches.}
        );
        return;
    }

    ## Check if the VEP api version and cache versions matches
    if ( any { not /$vep_version/xms } @vep_cache_version_folders ) {
        $log->fatal( q{Differing versions between '--vep_directory_path':}
              . $SPACE
              . $vep_directory_path
              . $SPACE
              . q{and '--vep_directory_cache':}
              . $SPACE
              . $vep_directory_cache );
        exit 1;
    }
    return 1;
}

sub _check_parameter_mandatory_keys_exits {

## Function : Check that mandatory keys exists
## Returns  :
## Arguments: $file_path          => Path to yaml file
##          : $mandatory_key_href => Hash with mandatory key {REF}
##          : $parameter_href     => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $mandatory_key_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        mandatory_key_href => {
            default     => {},
            required    => 1,
            store       => \$mandatory_key_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter ( keys %{$parameter_href} ) {

      MANDATORY_KEY:
        foreach my $mandatory_key ( keys %{$mandatory_key_href} ) {

            ## Mandatory key exists
            if ( not exists $parameter_href->{$parameter}{$mandatory_key} ) {

                say {*STDERR} q{Missing mandatory key: '}
                  . $mandatory_key
                  . q{' for parameter: '}
                  . $parameter
                  . q{' in file: '}
                  . $file_path
                  . $SINGLE_QUOTE
                  . $NEWLINE;
                croak();
            }
        }
    }
    return;
}

sub _check_parameter_keys {

## Function : Evaluate parameter keys in hash
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key_href       => Hash with mandatory key {REF}
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $key_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        key_href => {
            default     => {},
            required    => 1,
            store       => \$key_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter ( keys %{$parameter_href} ) {

      KEY:
        foreach my $key ( keys %{$key_href} ) {

            ## Key exists in parameter
            if ( exists $parameter_href->{$parameter}{$key} ) {

                ## Check key data type
                _check_parameter_data_type(
                    {
                        file_path      => $file_path,
                        key            => $key,
                        key_href       => $key_href,
                        parameter      => $parameter,
                        parameter_href => $parameter_href,
                    }
                );

                ## Evaluate key values
                _check_parameter_values(
                    {
                        file_path      => $file_path,
                        key            => $key,
                        key_href       => $key_href,
                        parameter      => $parameter,
                        parameter_href => $parameter_href,
                    }
                );
            }
        }
    }
    return;
}

sub _check_parameter_values {

## Function : Evaluate parameter key values
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key            => Hash with non  key
##          : $key_href       => Hash with  key {REF}
##          : $parameter      => Parameter
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $key;
    my $key_href;
    my $parameter;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        key =>
          { defined => 1, required => 1, store => \$key, strict_type => 1, },
        key_href => {
            default     => {},
            required    => 1,
            store       => \$key_href,
            strict_type => 1,
        },
        parameter => {
            defined     => 1,
            required    => 1,
            store       => \$parameter,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check value(s)
    if ( $key_href->{$key}{values} ) {

        my $value = $parameter_href->{$parameter}{$key};

        if ( not( any { $_ eq $value } @{ $key_href->{$key}{values} } ) ) {

            say {*STDERR} q{Found illegal value '}
              . $value
              . q{' for parameter: '}
              . $parameter
              . q{' in key: '}
              . $key
              . q{' in file: '}
              . $file_path
              . $SINGLE_QUOTE
              . $NEWLINE
              . q{Allowed entries: '}
              . join( q{', '}, @{ $key_href->{$key}{values} } )
              . $SINGLE_QUOTE
              . $NEWLINE;
            croak();
        }
    }
    return;
}

sub _check_parameter_data_type {

## Function : Check key data type
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key            => Hash with non key
##          : $key_href       => Hash with key {REF}
##          : $parameter      => Parameter
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $key;
    my $key_href;
    my $parameter;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        key =>
          { defined => 1, required => 1, store => \$key, strict_type => 1, },
        key_href => {
            default     => {},
            required    => 1,
            store       => \$key_href,
            strict_type => 1,
        },
        parameter => {
            defined     => 1,
            required    => 1,
            store       => \$parameter,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check data_type
    my $data_type = ref( $parameter_href->{$parameter}{$key} );

    ## Array or hash
    if ($data_type) {

        ## Wrong data_type
        if ( not $data_type eq $key_href->{$key}{key_data_type} ) {

            say {*STDERR} q{Found '}
              . $data_type
              . q{' but expected datatype '}
              . $key_href->{$key}{key_data_type}
              . q{' for parameter: '}
              . $parameter
              . q{' in key: '}
              . $key
              . q{' in file: '}
              . $file_path
              . $SINGLE_QUOTE
              . $NEWLINE;
            croak();
        }
    }
    elsif ( $key_href->{$key}{key_data_type} ne q{SCALAR} ) {

        ## Wrong data_type
        say {*STDERR} q{Found 'SCALAR' but expected datatype '}
          . $key_href->{$key}{key_data_type}
          . q{' for parameter: '}
          . $parameter
          . q{' in key: '}
          . $key
          . q{' in file: '}
          . $file_path
          . $SINGLE_QUOTE
          . $NEWLINE;
        croak();
    }
    return;
}

1;
