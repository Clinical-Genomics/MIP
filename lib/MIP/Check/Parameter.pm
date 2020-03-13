package MIP::Check::Parameter;

use 5.026;
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
use List::MoreUtils qw { all any uniq };

## MIPs lib/
use MIP::Constants qw{ $COMMA $DOLLAR_SIGN $DOT $LOG_NAME $NEWLINE $SINGLE_QUOTE $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.34;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_active_installation_parameters
      check_allowed_array_values
      check_infile_contain_sample_id
      check_infiles
      check_nist_file_exists
      check_nist_file_name
      check_nist_nist_id
      check_nist_sample_id
      check_nist_version
      check_prioritize_variant_callers
      check_recipe_fastq_compatibility
      check_sample_id_in_hash_parameter
      check_sample_id_in_hash_parameter_path
      check_select_file_contigs
      check_vep_custom_annotation
      check_vep_api_cache_versions
      check_vep_plugin
    };
}

sub check_active_installation_parameters {

## Function : Some active_parameter checks that are common to both installations. Returns "1" if all is OK
## Returns  : 1 or exit
## Arguments: $project_id => Project id
##          : sbatch_mode => Sbatch mode boolean

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $project_id;
    my $sbatch_mode;

    my $tmpl = {
        project_id => {
            store       => \$project_id,
            strict_type => 1,
        },
        sbatch_mode => {
            store       => \$sbatch_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check that a project id has been set if SBATCH mode
    if ( $sbatch_mode
        and not $project_id )
    {
        $log->fatal(
q{The parameter "project_id" must be set when a sbatch installation has been requested}
        );
        exit 1;
    }
    return 1;
}

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

sub check_nist_file_exists {

## Function : Check nist file path exists
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $log                   => Log object
##          : $nist_parameters_ref   => Nist parameters to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $nist_parameters_ref;

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
        nist_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$nist_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw { check_filesystem_objects_and_index_existance };

  NIST_PARAMETER:
    foreach my $nist_parameter ( @{$nist_parameters_ref} ) {

        # Alias
        my $nist_href = \%{ $active_parameter_href->{$nist_parameter} };

      NIST_VERSION:
        foreach my $nist_version ( keys %{$nist_href} ) {

          NIST_FILE:
            while ( my ( $nist_id, $file_path ) = each %{ $nist_href->{$nist_version} } )
            {
                ## Check path object exists
                check_filesystem_objects_and_index_existance(
                    {
                        object_name    => ( join q{=>}, ( $nist_version, $nist_id ) ),
                        object_type    => q{file},
                        parameter_name => $nist_parameter,
                        path           => $file_path,
                    }
                );

            }
        }
    }
    return 1;
}

sub check_nist_file_name {

## Function : Check nist file name is defined in nist parameters
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $log                   => Log object
##          : $nist_parameters_ref   => Nist parameters to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $nist_parameters_ref;

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
        nist_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$nist_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  NIST_PARAMETER:
    foreach my $nist_parameter ( @{$nist_parameters_ref} ) {

        # Alias
        my $nist_href = \%{ $active_parameter_href->{$nist_parameter} };

      NIST_VERSION:
        foreach my $nist_version ( keys %{$nist_href} ) {

          NIST_FILE:
            while ( my ( $nist_id, $file_name ) = each %{ $nist_href->{$nist_version} } )
            {

                ## Require that a file name is defined
                next NIST_FILE if ( defined $file_name );

                $log->fatal(
                    q{Please supply a file name for option: } . join q{=>},
                    ( $nist_parameter, $nist_version, $nist_id )
                );
                exit 1;

            }
        }
    }
    return 1;
}

sub check_nist_nist_id {

## Function : Check nist_ids contain supplied nist_id in nist parameters
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $log                   => Log object
##          : $nist_id_href          => Nist ids
##          : $nist_parameters_ref   => Nist parameters to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $nist_id_href;
    my $nist_parameters_ref;

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
        nist_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$nist_id_href,
            strict_type => 1,
        },
        nist_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$nist_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %seen;

  NIST_ID:
    while ( my ( $sample_id, $nist_id ) = each %{$nist_id_href} ) {

      NIST_PARAMETER:
        foreach my $nist_parameter ( @{$nist_parameters_ref} ) {

            # Alias
            my $nist_href = \%{ $active_parameter_href->{$nist_parameter} };

          NIST_VERSION:
            foreach my $nist_version ( keys %{$nist_href} ) {

                next NIST_VERSION if ( not exists $nist_href->{$nist_version}{$nist_id} );
                $seen{$nist_id}++;
            }
        }
        next NIST_ID if ( $seen{$nist_id} );

        $log->fatal(
            q{Supplied nist id: }
              . $nist_id
              . q{ for option --nist_id is not a defined nist_id supplied nist options: }
              . join $SPACE,
            @{$nist_parameters_ref}
        );
        exit 1;
    }
    return 1;
}

sub check_nist_sample_id {

## Function : Check nist_ids contain supplied sample_ids
## Returns  : 1
## Arguments: $log            => Log object
##          : $nist_id_href   => Nist ids
##          : $sample_ids_ref => Sample ids

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $nist_id_href;
    my $sample_ids_ref;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        nist_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$nist_id_href,
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

  SAMPLE_ID:
    while ( my ( $sample_id, $nist_id ) = each %{$nist_id_href} ) {

        ## Find supplied sample id
        next SAMPLE_ID if ( any { $_ eq $sample_id } @{$sample_ids_ref} );

        $log->fatal(
            q{Supplied sample id for option --nist_id }
              . q{ is not a supplied sample id: }
              . join $SPACE,
            @{$sample_ids_ref}
        );
        exit 1;
    }
    return 1;
}

sub check_nist_version {

## Function : Check nist_versions contain supplied nist_version in nist parameters
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $log                   => Log object
##          : $nist_parameters_ref   => Nist parameters to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $nist_parameters_ref;

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
        nist_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$nist_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Unpack
    my @nist_versions = @{ $active_parameter_href->{nist_versions} };

  NIST_PARAMETER:
    foreach my $nist_parameter ( @{$nist_parameters_ref} ) {

        # Alias
        my $nist_href = \%{ $active_parameter_href->{$nist_parameter} };

        ## Check that version exists in nist hashes
        next NIST_PARAMETER if ( all { exists $nist_href->{$_} } @nist_versions );

        $log->fatal(
            q{One or more nist versions }
              . ( join $SPACE, @nist_versions )
              . q{ does not exist in: }
              . join $SPACE,
            @{$nist_parameters_ref}
        );
        exit 1;
    }
    return 1;
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

        ## Only use first part of name
        my ($variant_caller_alias) = split /_/sxm, $variant_caller;
        push @variant_caller_aliases, $variant_caller_alias;

        ## Only active recipes
        if ( $active_parameter_href->{$variant_caller} ) {

            ## If variant caller alias is not part of priority order names
            if ( not any { $_ eq $variant_caller_alias } @priority_order_names ) {

                $log->fatal( $parameter_name
                      . q{ does not contain active variant caller: '}
                      . $variant_caller_alias
                      . $SINGLE_QUOTE );
                exit 1;
            }
        }
        else {
            ## Only NOT active recipes

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
            my $sample_id_value = $active_parameter_href->{$parameter_name}{$sample_id};
            my $is_mandatory    = $parameter_href->{$parameter_name}{mandatory};

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
        foreach my $sample_id_str ( keys %{ $active_parameter_href->{$parameter_name} } )
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

sub check_select_file_contigs {

## Function : Check that select file contigs is a subset of primary contigs
## Returns  :
## Arguments: $contigs_ref             => Primary contigs of the human genome reference
##          : $log                     => Log object to write to
##          : $select_file_contigs_ref => Select file contigs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $log;
    my $select_file_contigs_ref;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        select_file_contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$select_file_contigs_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ array_minus };

    ## Check that select file contigs are a subset of primary contigs
    my @unique_select_contigs =
      array_minus( @{$select_file_contigs_ref}, @{$contigs_ref} );
    if (@unique_select_contigs) {

        $log->fatal( q{Option 'vcfparser_select_file' contig(s): } . join $SPACE,
            @unique_select_contigs );
        $log->fatal(
            q{Is not a subset of the human genome reference contigs: } . join $SPACE,
            @{$contigs_ref} );
        exit 1;
    }
    return 1;
}

sub check_vep_custom_annotation {

## Function : Check VEP custom annotations options
## Returns  :
## Arguments: $log                 => Log object
##          : $vep_custom_ann_href => VEP custom annotation {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $vep_custom_ann_href;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        vep_custom_ann_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vep_custom_ann_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw { check_filesystem_objects_and_index_existance };

    ## Nothing to check
    return 0 if ( not keys %{$vep_custom_ann_href} );

  ANN:
    while ( my ( $ann, $value_href ) = each %{$vep_custom_ann_href} ) {

        my $err_msg = $ann . q{ Is not a hash ref for vep_custom_annotation};
        croak($err_msg) if ( ref $value_href ne q{HASH} );

        ## Check the VEP custom annotations options and that they have allowed values
        _check_vep_custom_annotation_options(
            {
                annotation             => $ann,
                custom_ann_option_href => $value_href,
                log                    => $log,
            }
        );

        ## Check path object exists
        check_filesystem_objects_and_index_existance(
            {
                object_name    => $ann,
                object_type    => q{file},
                parameter_name => q{vep_custom_annotation},
                path           => $value_href->{path},
            }
        );
    }
    return 1;
}

sub check_vep_api_cache_versions {

## Function : Compare VEP API and VEP chache versions. Exit if non-match
## Returns  :
## Arguments: $vep_binary_path     => VEP binary path {REF}
##          : $vep_directory_cache => VEP cache directory path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_directory_cache;

    ## Default(s)
    my $vep_binary_path;

    my $tmpl = {
        vep_binary_path => {
            default     => q{vep},
            defined     => 1,
            store       => \$vep_binary_path,
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
    use MIP::Get::Parameter qw{ get_vep_version };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $vep_version = get_vep_version(
        {
            vep_bin_path => $vep_binary_path,
        }
    );

    ## Check that a version number was picked up
    if ( not $vep_version ) {
        $log->warn(
q{Could not retrieve VEP version. Skipping checking that VEP api and cache matches.}
        );
        return;
    }

    ## Get absolute path to VEP species cache as it is commonly linked
    my $vep_cache_species_dir_path =
      _get_vep_cache_species_dir_path( { vep_directory_cache => $vep_directory_cache, } );

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
    my @vep_cache_version_folders = $rule->in($vep_cache_species_dir_path);

    ## Check that
    if ( not @vep_cache_version_folders ) {
        $log->warn(
q{Could not retrieve VEP cache version. Skipping checking that VEP api and cache matches.}
        );
        return;
    }

    ## Check if the VEP api version and cache versions matches
    if ( not any { /$vep_version/xms } @vep_cache_version_folders ) {

        $log->fatal( q{Differing versions between 'VEP API':}
              . $SPACE
              . $vep_version
              . $SPACE
              . q{and '--vep_directory_cache':}
              . $SPACE
              . $vep_directory_cache );
        exit 1;
    }
    return 1;
}

sub check_vep_plugin {

## Function : Check VEP plugin options
## Returns  : 0 or 1
## Arguments: $log             => Log object
##          : $parameter_name  => Parameter name
##          : $vep_plugin_href => VEP plugin annotation {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $parameter_name;
    my $vep_plugin_href;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_name => {
            defined  => 1,
            required => 1,
            store    => \$parameter_name,
        },
        vep_plugin_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vep_plugin_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw { check_filesystem_objects_and_index_existance };

    ## Nothing to check
    return 0 if ( not keys %{$vep_plugin_href} );

  PLUGIN:
    while ( my ( $plugin, $value_href ) = each %{$vep_plugin_href} ) {

        my $err_msg = $plugin . q{ Is not a hash ref for vep_plugin};
        croak($err_msg) if ( ref $value_href ne q{HASH} );

        next PLUGIN if ( not exists $value_href->{path} );

        next PLUGIN if ( not exists $value_href->{exists_check} );

        ## Check path object exists
        check_filesystem_objects_and_index_existance(
            {
                object_name    => $plugin,
                object_type    => $value_href->{exists_check},
                parameter_name => $parameter_name,
                path           => $value_href->{path},
            }
        );
    }
    return 1;
}

sub check_recipe_fastq_compatibility {

## Function : Check that the recipe is compatible with the fastq sequence modes. Turn of downstream applications otherwise
## Returns  :
## Arguments: $active_parameter_href   => Active parameter hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name             => Recipe name
##          : $sample_info_href        => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten arguments
    my $active_parameter_href;
    my $infile_lane_prefix_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            required    => 1,
            defined     => 1,
            store       => \$recipe_name,
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

    use MIP::Dependency_tree
      qw{ get_recipe_dependency_tree_chain get_recipes_for_dependency_tree_chain };
    use MIP::Sample_info qw{ get_sequence_run_type };
    use MIP::Set::Parameter qw{ set_recipe_mode };

    ## Check if program is going to run
    return if ( $active_parameter_href->{$recipe_name} == 0 );

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $is_compatible = 1;

    ## Get sequence run modes
  SAMPLE_ID:
    foreach my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

        my %sequence_run_type = get_sequence_run_type(
            {
                infile_lane_prefix_href => $infile_lane_prefix_href,
                sample_id               => $sample_id,
                sample_info_href        => $sample_info_href,
            }
        );

        ## Turn of recipe if multiple sequence types are present
        if ( uniq( values %sequence_run_type ) > 1 ) {
            $is_compatible = 0;
        }
    }

    if ( not $is_compatible ) {

        ## Turn of current recipe and downstream recipes
        $log->warn(q{Multiple sequence run types detected});
        $log->warn(qq{Turning off $recipe_name and downstream recipes});

        my $recipe_chain;
        get_recipe_dependency_tree_chain(
            {
                recipe               => $recipe_name,
                dependency_tree_href => $parameter_href->{dependency_tree_href},
                chain_id_ref         => \$recipe_chain,
            }
        );

        my @chain_recipes = get_recipes_for_dependency_tree_chain(
            {
                dependency_tree_href    => $parameter_href->{dependency_tree_href},
                chain_initiation_point  => $recipe_chain,
                recipe_initiation_point => $recipe_name,
            }
        );

        ## Turn of recipes
        set_recipe_mode(
            {
                active_parameter_href => $active_parameter_href,
                recipes_ref           => \@chain_recipes,
                mode                  => 0,
            }
        );
    }

    return $is_compatible;
}

sub _check_vep_custom_annotation_options {

## Function : Check the VEP custom annotations options are defined and with allowed values
## Returns  :
## Arguments: $annotation             => Annotation
##          : $custom_ann_option_href => Custom annotation options
##          : $log                    => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotation;
    my $custom_ann_option_href;
    my $log;

    my $tmpl = {
        annotation => {
            defined     => 1,
            required    => 1,
            store       => \$annotation,
            strict_type => 1,
        },
        custom_ann_option_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$custom_ann_option_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %check_vep_annotations = (
        annotation_type          => { allow => [qw{ exact overlap }], },
        file_type                => { allow => [qw{ bed gff gtf vcf bigwig }], },
        force_report_coordinates => { allow => [ 0, 1 ], },
    );

    ## Check required keys
    my @required_options = (qw{ key });
  REQ_OPTION:
    foreach my $required_option (@required_options) {

        if (   not exists $custom_ann_option_href->{$required_option}
            or not defined $custom_ann_option_href->{$required_option} )
        {

            $log->fatal( q{Vep custom annotation option hash: }
                  . $annotation
                  . q{ lacks required option }
                  . $required_option );
            exit 1;
        }
    }

    ## Check allowed options for annotation
  OPTION:
    foreach my $option ( keys %{$custom_ann_option_href} ) {

        ## Allow anything defined
        next OPTION if ( $option eq q{key} );

        next OPTION if ( $option eq q{path} );

        next OPTION if ( $option eq q{vcf_fields} );

        next OPTION
          if (
            any { $_ eq $custom_ann_option_href->{$option} }
            @{ $check_vep_annotations{$option}{allow} }
          );

        $log->fatal( q{Vep custom annotation option hash: }
              . $annotation
              . q{ has a not allowed option value '}
              . $option . q{ => }
              . $custom_ann_option_href->{$option} );
        $log->fatal(
            q{Allowed options are: } . join $SPACE,
            @{ $check_vep_annotations{$option}{allow} }
        );
        exit 1;
    }

    return 1;
}

sub _get_vep_cache_species_dir_path {

## Function : Get the vep cache species dir path
## Returns  : $vep_cache_species_dir_path
## Arguments: $vep_directory_cache => VEP cache directory path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_directory_cache;

    my $tmpl = {
        vep_directory_cache => {
            defined     => 1,
            required    => 1,
            store       => \$vep_directory_cache,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Cwd qw{ abs_path };

    my @vep_species_cache = qw{ homo_sapiens homo_sapiens_merged };
    my $vep_cache_species_dir_path;

  SPECIES_CACHE:
    foreach my $species_cache (@vep_species_cache) {

        ## Get species specific cache dir
        $vep_cache_species_dir_path =
          abs_path( catdir( $vep_directory_cache, $species_cache ) );
        last if ( -e $vep_cache_species_dir_path );
    }
    return $vep_cache_species_dir_path;
}

1;
