package MIP::Pedigree;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ make_path };
use List::Util qw{ none };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error};
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use List::MoreUtils qw { any };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $LOG_NAME $NEWLINE $SPACE $TAB $UNDERSCORE };

BEGIN {
    use base qw{ Exporter };
    require Exporter;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_founder_id
      check_pedigree_mandatory_key
      check_pedigree_sample_allowed_values
      check_pedigree_vs_user_input_sample_ids
      create_fam_file
      get_is_trio
      gatk_pedigree_flag
      has_duo
      has_trio
      is_sample_proband_in_trio
      parse_pedigree
      parse_yaml_pedigree_file
      set_active_parameter_pedigree_keys
      set_pedigree_capture_kit_info
      set_pedigree_case_info
      set_pedigree_phenotype_info
      set_pedigree_sample_info
    };
}

## Constants
Readonly my $DUO_MEMBERS_COUNT  => 2;
Readonly my $TRIO_MEMBERS_COUNT => 3;

sub check_founder_id {

## Function : Check that founder_ids are included in the pedigree info
## Returns  :
## Arguments: $active_sample_ids_ref => Array of pedigree samples {REF}
##          : $pedigree_href         => Pedigree info {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_sample_ids_ref;
    my $pedigree_href;

    my $tmpl = {
        active_sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$active_sample_ids_ref,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

  SAMPLE:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        my @founders =
          ( $pedigree_sample_href->{father}, $pedigree_sample_href->{mother} );

      FOUNDER:
        foreach my $founder (@founders) {

            ## No founder info
            next FOUNDER if ( not $founder );

            ## If founder is part of analysis
            next FOUNDER
              if ( any { $_ eq $founder } @{$active_sample_ids_ref} );

            $log->fatal( q{Could not find founder sample_id: }
                  . $founder
                  . q{ from pedigree file in current analysis} );
            exit 1;
        }
    }
    return 1;
}

sub check_pedigree_mandatory_key {

## Function : Check that the pedigree case mandatory keys are present
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $pedigree_file_path    => Pedigree file path
##          : $pedigree_href         => YAML pedigree info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $pedigree_file_path;
    my $pedigree_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        pedigree_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$pedigree_file_path,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @mandatory_case_keys   = qw{ case samples };
    my @mandatory_sample_keys = qw{ sample_id father mother sex phenotype };

    ## Add analysis specific mandatory keys
    if ( $active_parameter_href->{dna_vcf_file} ) {

        push @mandatory_sample_keys, q{dna_sample_id};
    }

    ### Case
    ## Check mandatory case keys
  MANDATORY_KEY:
    foreach my $key (@mandatory_case_keys) {

        next MANDATORY_KEY if ( exists $pedigree_href->{$key} );

        $log->fatal( q{Pedigree file: }
              . $pedigree_file_path
              . q{ cannot find mandatory key: }
              . $key
              . q{ in file} );
        exit 1;
    }

    ### Sample
    ## Check mandatory sample keys
  SAMPLE_KEY:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        ## Check that we find mandatory sample keys
      MANDATORY_KEY:
        foreach my $key (@mandatory_sample_keys) {

            next MANDATORY_KEY if ( exists $pedigree_sample_href->{$key} );

            $log->fatal( q{Pedigree file: }
                  . $pedigree_file_path
                  . q{ cannot find mandatory sample key: }
                  . $key
                  . q{ in file} );
            exit 1;
        }
    }
    return 1;
}

sub check_pedigree_sample_allowed_values {

## Function : Check that the pedigree sample key values are allowed
## Returns  :
## Arguments: $pedigree_file_path => Pedigree file path
##          : $pedigree_href      => YAML pedigree info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $pedigree_file_path;
    my $pedigree_href;

    my $tmpl = {
        pedigree_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$pedigree_file_path,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %allowed_values = (
        analysis_type => [qw{ dragen_rd_dna panel vrn wes wgs wts }],
        phenotype     => [qw{ affected unaffected unknown }],
        sex           => [qw{ female male other unknown }],
    );

    ## Check input to sample_info hash at sample level
  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

      SAMPLE_KEY:
        foreach my $key ( keys %{$pedigree_sample_href} ) {

            ## No defined allowed values
            next SAMPLE_KEY if ( not exists $allowed_values{$key} );

            ## If element is not part of array
            next SAMPLE_KEY
              if ( any { $_ eq $pedigree_sample_href->{$key} } @{ $allowed_values{$key} } );

            $log->fatal( q{Pedigree file: } . $pedigree_file_path );
            $log->fatal( q{For key: }
                  . $key
                  . q{ found illegal value: '}
                  . $pedigree_sample_href->{$key}
                  . q{' allowed values are '}
                  . join( q{' '}, @{ $allowed_values{$key} } )
                  . q{'} );
            $log->fatal(q{Please correct the entry before analysis});
            $log->fatal(q{Aborting run});
            exit 1;
        }
    }
    return 1;
}

sub check_pedigree_vs_user_input_sample_ids {

## Function : Check pedigree sample ids and user input sample ids match
## Returns  :
## Arguments: $pedigree_file_path        => Pedigree file path
##          : $pedigree_sample_ids_ref   => Pedigree sample ids
##          : $user_input_sample_ids_ref => User input sample ids

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $pedigree_file_path;
    my $pedigree_sample_ids_ref;
    my $user_input_sample_ids_ref;

    my $tmpl = {
        pedigree_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$pedigree_file_path,
            strict_type => 1,
        },
        pedigree_sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$pedigree_sample_ids_ref,
            strict_type => 1,
        },
        user_input_sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$user_input_sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

  SAMPLE_ID:
    foreach my $sample_id ( @{$user_input_sample_ids_ref} ) {

        ## If element is not part of array
        if ( not any { $_ eq $sample_id } @{$pedigree_sample_ids_ref} ) {

            $log->fatal( q{File: }
                  . $pedigree_file_path
                  . q{ provided sample_id: }
                  . $sample_id
                  . q{ is not present in pedigree file} );
            exit 1;
        }
    }
    return 1;
}

sub create_fam_file {

## Function : Create ".fam" file to be used in variant calling analyses
## Returns  :
## Arguments: $case_id          => Case_id
##          : $execution_mode   => Either system (direct) or via sbatch
##          : $fam_file_path    => Case file path
##          : $filehandle       => Filehandle to write to {Optional unless execution_mode=sbatch}
##          : $include_header   => Include header ("1") or not ("0")
##          : $sample_ids_ref   => Sample ids {REF}
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $fam_file_path;
    my $filehandle;
    my $sample_ids_ref;
    my $sample_info_href;

    ## Default(s)
    my $execution_mode;
    my $include_header;

    my $tmpl = {
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        execution_mode => {
            allow       => [qw{ sbatch system }],
            default     => q{sbatch},
            store       => \$execution_mode,
            strict_type => 1,
        },
        fam_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$fam_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        include_header => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$include_header,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
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

    use MIP::Sample_info qw{ set_in_sample_info };

    my @fam_lines = _build_fam_file_line(
        {
            case_id          => $case_id,
            include_header   => $include_header,
            sample_ids_ref   => $sample_ids_ref,
            sample_info_href => $sample_info_href,
        }
    );

    _write_fam_file(
        {
            execution_mode => $execution_mode,
            fam_file_path  => $fam_file_path,
            fam_lines_ref  => \@fam_lines,
            filehandle     => $filehandle,
        }
    );

    return;
}

sub get_is_trio {

## Function  : Detect case constellation based on pedigree file
## Returns   : undef | 1
## Arguments : $active_parameter_href => Active parameters for this analysis hash {REF}
##           : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
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

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    if ( scalar @{ $active_parameter_href->{sample_ids} } == 1 ) {

        $log->info( q{Found single sample: } . $active_parameter_href->{sample_ids}[0] );
        return;
    }
    elsif ( scalar @{ $active_parameter_href->{sample_ids} } == $TRIO_MEMBERS_COUNT ) {

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Unpack
            my $father_id =
              $sample_info_href->{sample}{$sample_id}{father};
            my $mother_id =
              $sample_info_href->{sample}{$sample_id}{mother};
            my $is_trio = _parse_trio_members(
                {
                    father_id      => $father_id,
                    mother_id      => $mother_id,
                    sample_id      => $sample_id,
                    sample_ids_ref => $active_parameter_href->{sample_ids},
                }
            );
            ## Return if a trio is found
            return $is_trio if ($is_trio);
        }
    }
    return;
}

sub has_duo {

## Function  : Check if case has a parent-child duo and child is affected
## Returns   : 0 | 1
## Arguments : $active_parameter_href => Active parameters for this analysis hash {REF}
##           : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
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

    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    ## Has two samples
    return 0
      if ( scalar @{ $active_parameter_href->{sample_ids} } != $DUO_MEMBERS_COUNT );

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        my %sample_attributes = get_pedigree_sample_id_attributes(
            {
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        ## Find a child
        next SAMPLE_ID if ( not( $sample_attributes{father} or $sample_attributes{mother} ) );

        return 1 if ( $sample_attributes{phenotype} eq q{affected} );
    }
    return 0;
}

sub has_trio {

## Function  : Check if case has trio
## Returns   : 0 | 1
## Arguments : $active_parameter_href => Active parameters for this analysis hash {REF}
##           : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
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

    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    ## At least three samples
    return 0
      if ( scalar @{ $active_parameter_href->{sample_ids} } < $TRIO_MEMBERS_COUNT );

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        my $mother = get_pedigree_sample_id_attributes(
            {
                attribute        => q{mother},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        my $father = get_pedigree_sample_id_attributes(
            {
                attribute        => q{father},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        ## Find a child
        next SAMPLE_ID if ( not $father or not $mother );

        my $phenotype = get_pedigree_sample_id_attributes(
            {
                attribute        => q{phenotype},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        return 1 if ( $phenotype eq q{affected} );
    }
    return 0;
}

sub gatk_pedigree_flag {

## Function : Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
## Returns  : %command
## Arguments: $fam_file_path            => The case file path
##          : $pedigree_validation_type => The pedigree validation strictness level

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fam_file_path;

    ## Default(s)
    my $pedigree_validation_type;

    my $tmpl = {
        fam_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$fam_file_path,
            strict_type => 1,
        },
        pedigree_validation_type => {
            allow       => [qw{ SILENT STRICT }],
            default     => q{SILENT},
            store       => \$pedigree_validation_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $parent_counter;
    my $pq_parent_counter = _build_parent_child_counter_regexp(
        {
            case_member => q{parent},
        }
    );

    my $child_counter;
    my $pq_child_counter = _build_parent_child_counter_regexp(
        {
            case_member => q{child},
        }
    );

    my %command;

    ## Count the number of parents
    $parent_counter = `$pq_parent_counter $fam_file_path`;

    ## Count the number of children
    $child_counter = `$pq_child_counter $fam_file_path`;

    if ( $parent_counter > 0 ) {

        $command{pedigree_validation_type} = $pedigree_validation_type;

        ## Pedigree files for samples
        $command{pedigree} = $fam_file_path;
    }
    return %command;
}

sub is_sample_proband_in_trio {

## Function : Check if sample id has an affected or unknown phenotype and is child in trio
## Returns  : 0 | 1
## Arguments: $only_affected    => Control whether the proband needs to marked as affected
##          : $sample_id        => Sample id
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $only_affected;

    my $tmpl = {
        only_affected => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$only_affected,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    use MIP::Sample_info qw{ get_family_member_id get_pedigree_sample_id_attributes };

    ## There has to be a trio
    return 0 if ( not $sample_info_href->{has_trio} );

    ## Get phenotype
    my $phenotype = get_pedigree_sample_id_attributes(
        {
            attribute        => q{phenotype},
            sample_id        => $sample_id,
            sample_info_href => $sample_info_href,
        }
    );

    ## Sample_id needs to be affected
    return 0 if ( $only_affected && ( $phenotype eq q{unaffected} ) );

    ## Get family hash
    my %family_member_id = get_family_member_id( { sample_info_href => $sample_info_href } );

    ## Check if the sample is an affected child
    return 0 if ( none { $_ eq $sample_id } @{ $family_member_id{children} } );

    return 1;
}

sub parse_pedigree {

## Function : Parse pedigree info pedigree file and updates pedigree parameters
##            with user supplied info
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $pedigree_file_path    => Pedigree file path from cmd or sample_info
##          : $parameter_href        => Parameter hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $pedigree_file_path;
    my $parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        pedigree_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$pedigree_file_path,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
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

    use MIP::Io::Read qw{ read_from_file };

    return if ( not defined $pedigree_file_path );

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Load parameters for this run or if not provided on cmd - load
    ## previous run from sample_info_file
    my %pedigree = read_from_file(
        {
            format => q{yaml},
            path   => $pedigree_file_path,
        }
    );

    $log->info( q{Loaded: } . $pedigree_file_path );

    ## Check pedigree data for allowed entries and correct format
    ## Add data to sample_info depending on user info
    parse_yaml_pedigree_file(
        {
            active_parameter_href => $active_parameter_href,
            pedigree_file_path    => $pedigree_file_path,
            parameter_href        => $parameter_href,
            pedigree_href         => \%pedigree,
            sample_info_href      => $sample_info_href,
        }
    );
    return 1;
}

sub parse_yaml_pedigree_file {

## Function : Parse case info in YAML pedigree file. Check pedigree data for
##            allowed entries and correct format. Add data to sample_info and
##            active_parameter depending on user info.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $pedigree_file_path    => Pedigree file path
##          : $parameter_href        => Parameter hash {REF}
##          : $pedigree_href         => Pedigree hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $pedigree_file_path;
    my $parameter_href;
    my $pedigree_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        pedigree_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$pedigree_file_path,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
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

    use MIP::Active_parameter qw{ get_user_supplied_pedigree_parameter };
    use MIP::Parameter qw{ get_capture_kit };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Use to collect which sample_ids have used a certain capture_kit
    my $case_id = $pedigree_href->{case};

    ## User supplied sample_ids via cmd or config
    my @user_input_sample_ids;

    ## Check pedigree mandatory keys
    check_pedigree_mandatory_key(
        {
            active_parameter_href => $active_parameter_href,
            pedigree_file_path    => $pedigree_file_path,
            pedigree_href         => $pedigree_href,
        }
    );

    ## Check that supplied YAML pedigree case_id and cmd match
    if ( $pedigree_href->{case} ne $active_parameter_href->{case_id} ) {

        $log->fatal( q{Pedigree file: }
              . $pedigree_file_path
              . q{ for  pedigree case_id: '}
              . $pedigree_href->{case}
              . q{' and supplied case: '}
              . $active_parameter_href->{case_id}
              . q{' does not match} );
        exit 1;
    }

    set_pedigree_case_info(
        {
            pedigree_href    => $pedigree_href,
            sample_info_href => $sample_info_href,
        }
    );

    ### Check sample keys values
    check_pedigree_sample_allowed_values(
        {
            pedigree_file_path => $pedigree_file_path,
            pedigree_href      => $pedigree_href,
        }
    );

    ## Get potential input from user
    my %is_user_supplied = get_user_supplied_pedigree_parameter(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    if ( $is_user_supplied{sample_ids}
        and exists $active_parameter_href->{sample_ids} )
    {

        ## Set cmd or config supplied sample_ids
        @user_input_sample_ids = @{ $active_parameter_href->{sample_ids} };
    }

    my @pedigree_sample_ids = set_pedigree_sample_info(
        {
            active_parameter_href     => $active_parameter_href,
            is_user_supplied_href     => \%is_user_supplied,
            pedigree_href             => $pedigree_href,
            sample_info_href          => $sample_info_href,
            user_input_sample_ids_ref => \@user_input_sample_ids,
        }
    );

    if ( not $is_user_supplied{sample_ids} ) {

        ## Lexiographical sort to determine the correct order of ids indata
        @{ $active_parameter_href->{sample_ids} } =
          sort @{ $active_parameter_href->{sample_ids} };
    }
    else {

        ## Check that CLI supplied sample_id exists in pedigree
        check_pedigree_vs_user_input_sample_ids(
            {
                pedigree_file_path        => $pedigree_file_path,
                pedigree_sample_ids_ref   => \@pedigree_sample_ids,
                user_input_sample_ids_ref => \@user_input_sample_ids,
            }
        );
    }

    ## Add phenotype to dynamic parameters
    set_pedigree_phenotype_info(
        {
            parameter_href => $parameter_href,
            pedigree_href  => $pedigree_href,
        }
    );

    set_active_parameter_pedigree_keys(
        {
            active_parameter_href => $active_parameter_href,
            is_user_supplied_href => \%is_user_supplied,
            pedigree_href         => $pedigree_href,
            sample_info_href      => $sample_info_href,
        }
    );

    set_pedigree_capture_kit_info(
        {
            active_parameter_href => $active_parameter_href,
            is_user_supplied_href => \%is_user_supplied,
            parameter_href        => $parameter_href,
            pedigree_href         => $pedigree_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Check that founder_ids are included in the pedigree info and the analysis run
    check_founder_id(
        {
            active_sample_ids_ref => \@{ $active_parameter_href->{sample_ids} },
            pedigree_href         => $pedigree_href,
        }
    );

    return 1;
}

sub set_active_parameter_pedigree_keys {

## Function : Set the pedigree case keys and values
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $is_user_supplied_href => User supplied parameter switch {REF}
##          : $pedigree_href         => YAML pedigree info hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $is_user_supplied_href;
    my $pedigree_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        is_user_supplied_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$is_user_supplied_href,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
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

    use MIP::Active_parameter qw{ set_pedigree_sample_id_parameter };
    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    my @pedigree_keys = qw{ analysis_type expected_coverage time_point };

  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        ## Unpack
        my $sample_id = $pedigree_sample_href->{sample_id};

      PEDIGREE_KEY:
        foreach my $pedigree_key (@pedigree_keys) {

            my $pedigree_value = get_pedigree_sample_id_attributes(
                {
                    attribute        => $pedigree_key,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            ## Key was not set in pedigree (this run or previous run)
            ## Hence not stored in sample_info
            next PEDIGREE_KEY if ( not defined $pedigree_value );

            ## User supplied info on cmd or config
            next PEDIGREE_KEY if ( $is_user_supplied_href->{$pedigree_key} );

            ## Add value for sample_id using pedigree info
            set_pedigree_sample_id_parameter(
                {
                    active_parameter_href => $active_parameter_href,
                    pedigree_key          => $pedigree_key,
                    pedigree_value        => $pedigree_value,
                    sample_id             => $sample_id,
                }
            );
        }
    }
    return;
}

sub set_pedigree_capture_kit_info {

## Function : Add and set capture kit for each individual
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $pedigree_href         => YAML pedigree info hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}
##          : $is_user_supplied_href => User supplied info switch {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $pedigree_href;
    my $sample_info_href;
    my $is_user_supplied_href;

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
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        is_user_supplied_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$is_user_supplied_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ set_exome_target_bed };
    use MIP::Parameter qw{ get_capture_kit };
    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    my %exom_target_bed_file_sample_id_map;

  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        ## Unpack
        my $sample_id   = $pedigree_sample_href->{sample_id};
        my $capture_kit = get_pedigree_sample_id_attributes(
            {
                attribute        => q{capture_kit},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        ## No recorded capture kit from pedigree or previous run
        next SAMPLE_HREF if ( not $capture_kit );

        ## Return a capture kit depending on user info
        my $exome_target_bed_file = get_capture_kit(
            {
                capture_kit                => $capture_kit,
                supported_capture_kit_href => $parameter_href->{supported_capture_kit}{default},
                is_set_by_user             => $is_user_supplied_href->{exome_target_bed},
            }
        );

        ## No capture kit returned
        next SAMPLE_HREF if ( not $exome_target_bed_file );

        push @{ $exom_target_bed_file_sample_id_map{$exome_target_bed_file} }, $sample_id;
    }

  BED_FILE:
    foreach my $exome_target_bed_file ( keys %exom_target_bed_file_sample_id_map ) {

        ## We have read capture kits from pedigree and
        ## need to transfer to active_parameters
        my $sample_id_string = join $COMMA,
          @{ $exom_target_bed_file_sample_id_map{$exome_target_bed_file} };

        set_exome_target_bed(
            {
                active_parameter_href => $active_parameter_href,
                exome_target_bed_file => $exome_target_bed_file,
                sample_id_string      => $sample_id_string,
            }
        );
    }
    return;
}

sub set_pedigree_case_info {

## Function : Get the pedigree case keys and values
## Returns  :
## Arguments: $pedigree_href    => Pedigree info hash {REF}
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $pedigree_href;
    my $sample_info_href;

    my $tmpl = {
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
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

    ### Add key and values for case level info
  KEY:
    foreach my $key ( keys %{$pedigree_href} ) {

        ## Do not add sample level info
        next KEY if ( $key eq q{samples} );

        $sample_info_href->{$key} = $pedigree_href->{$key};
    }
    return;
}

sub set_pedigree_sample_info {

## Function : Set the pedigree sample keys and values
## Returns  :
## Arguments: $active_parameter_href     => Active parameters for this analysis hash {REF}
##          : $pedigree_href             => YAML pedigree info hash {REF}
##          : $sample_info_href          => Info on samples and case hash {REF}
##          : $is_user_supplied_href     => User supplied parameter switch {REF}
##          : $user_input_sample_ids_ref => User supplied sample_ids via cmd or config

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $pedigree_href;
    my $sample_info_href;
    my $is_user_supplied_href;
    my $user_input_sample_ids_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        is_user_supplied_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$is_user_supplied_href,
            strict_type => 1,
        },
        user_input_sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$user_input_sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @pedigree_sample_ids;

    ### Add values sample level info
  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        # Unpack
        my $sample_id = $pedigree_sample_href->{sample_id};

        ## Save pedigree sample_id info for later check
        push @pedigree_sample_ids, $sample_id;

        ## No sample_ids supplied by user
        if ( not $is_user_supplied_href->{sample_ids} ) {

            ## Save sample_id info for analysis
            push @{ $active_parameter_href->{sample_ids} }, $sample_id;

            ## Add input to sample_info hash for at sample level
            _add_pedigree_sample_info(
                {
                    pedigree_sample_href => $pedigree_sample_href,
                    sample_id            => $sample_id,
                    sample_info_href     => $sample_info_href,
                }
            );
            next SAMPLE_HREF;
        }

        ### Sample_ids supplied by user
        ## Update sample_id with pedigree info for user supplied sample_ids
        if ( any { $_ eq $sample_id } @{$user_input_sample_ids_ref} ) {

            ## Set pedigree sample info
            _add_pedigree_sample_info(
                {
                    pedigree_sample_href => $pedigree_sample_href,
                    sample_id            => $sample_id,
                    sample_info_href     => $sample_info_href,
                }
            );
        }
    }
    return @pedigree_sample_ids;
}

sub set_pedigree_phenotype_info {

## Function : Store phenotype in dynamic parameters.
##          : Reformats pedigree phenotype to plink format before adding
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}
##          : $pedigree_href  => YAML pedigree info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $pedigree_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        # Unpack
        my $sample_id = $pedigree_sample_href->{sample_id};
        my $phenotype = $pedigree_sample_href->{phenotype};

        next SAMPLE_HREF if ( not defined $phenotype );

        push @{ $parameter_href->{cache}{$phenotype} }, $sample_id;
    }
    return;
}

sub _add_pedigree_sample_info {

## Function : Add pedigree sample level keys and values to sample info
## Returns  :
## Arguments: $pedigree_sample_href => Pedigree sample info hash {REF}
##          : $sample_id            => Sample ID
##          : $sample_info_href     => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $pedigree_sample_href;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        pedigree_sample_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_sample_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    ## Add input to sample_info hash at sample level
  PEDIGREE_SAMPLE_KEY:
    foreach my $key ( keys %{$pedigree_sample_href} ) {

        $sample_info_href->{sample}{$sample_id}{$key} =
          $pedigree_sample_href->{$key};
    }
    return;
}

sub _build_parent_child_counter_regexp {

## Function : Create regexp to count the number of parents / children
## Returns  : $regexp
## Arguments: $case_member => Parent or child

    my ($arg_href) = @_;

    ## Flatten argument
    my $case_member;

    my $tmpl = {
        case_member => {
            defined  => 1,
            required => 1,
            store    => \$case_member,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    my $regexp = q?perl -ne '?;

    ## Add parent or child, according to which one needs to be counted
    $regexp .= q?my $? . $case_member . $UNDERSCORE . q?counter=0; ?;

    ## Split line around tab if it's not a comment
    $regexp .=
q?while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] eq 0) || ($line[3] eq 0) ) ?;

    ## Increment the counter
    $regexp .= q?{ $? . $case_member . q?++} } } print $? . $case_member . q?; last;'?;

    return $regexp;
}

sub _parse_trio_members {

## Function  : Parse trio constellation
## Returns   : %trio
## Arguments : $father_id      => Potential father
##           : $mother_id      => Potential mother
##           : $sample_id      => Sample under investigation
##           : $sample_ids_ref => Sample_ids in current analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $father_id;
    my $mother_id;
    my $sample_id;
    my $sample_ids_ref;

    my $tmpl = {
        father_id => {
            required => 1,
            store    => \$father_id,
        },
        mother_id => {
            required => 1,
            store    => \$mother_id,
        },
        sample_id => {
            required => 1,
            store    => \$sample_id,
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

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %trio;

    ## Not a child
    return if ( not $father_id and not $mother_id );

    ## Sample_id must be a child
    $trio{child} = $sample_id;

    my %parent = (
        father => $father_id,
        mother => $mother_id,
    );

  PARENT:
    while ( my ( $parent_role, $parent_id ) = each %parent ) {

        ## If parent is present in current analysis
        if ( any { $_ eq $parent_id } @{$sample_ids_ref} ) {

            ## Set as parents
            $trio{$parent_role} = $parent_id;
        }
    }
    if ( scalar( keys %trio ) == $TRIO_MEMBERS_COUNT ) {

        $log->info( q{Found trio: Child = "}
              . $trio{child}
              . q{", Father = "}
              . $trio{father}
              . q{", Mother = "}
              . $trio{mother}
              . q{"}, );
        return 1;
    }
    return;
}

sub _build_fam_file_line {

## Function : Build fam file lines for each sample
## Returns  : @fam_lines
## Arguments: $case_id          => Case_id
##          : $include_header   => Include header ("1") or not ("0")
##          : $sample_ids_ref   => Sample ids {REF}
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $sample_ids_ref;
    my $sample_info_href;

    ## Default(s)
    my $include_header;

    my $tmpl = {
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        include_header => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$include_header,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
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

    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    my @fam_headers = ( q{#family_id}, qw{ sample_id father mother sex phenotype } );
    my @fam_lines;

    ## Add fam file headers
    if ($include_header) {

        push @fam_lines, join $TAB, @fam_headers;
    }

  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        my %pedigree_sample_attribute = get_pedigree_sample_id_attributes(
            {
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        ## Build sample line string
        my $sample_line = $case_id;

      HEADER:
        foreach my $header (@fam_headers) {

            next HEADER if ( not defined $pedigree_sample_attribute{$header} );

            my $pedigree_attribute = $pedigree_sample_attribute{$header};
            $pedigree_attribute = _convert_to_fam_format(
                {
                    header             => $header,
                    pedigree_attribute => $pedigree_attribute,
                }
            );

            $sample_line .= $TAB . $pedigree_attribute;
        }
        push @fam_lines, $sample_line;
    }
    return @fam_lines;
}

sub _convert_to_fam_format {

## Function : Convert pedigree attribute to fam format
## Returns  : fam format pedigree attribute or $pedigree_attribute
## Arguments: $header             => Pedigree header
##          : $pedigree_attribute => Sample pedigree attribute

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $header;
    my $pedigree_attribute;

    my $tmpl = {
        header => {
            defined     => 1,
            required    => 1,
            store       => \$header,
            strict_type => 1,
        },
        pedigree_attribute => {
            defined     => 1,
            required    => 1,
            store       => \$pedigree_attribute,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %fam_conversion = (
        phenotype => {
            affected   => 2,
            other      => q{other},
            unaffected => 1,
            unknown    => 0,
        },
        sex => {
            female  => 2,
            male    => 1,
            other   => q{other},
            unknown => q{other},
        },
    );

    return ${fam_conversion}{$header}{$pedigree_attribute} // $pedigree_attribute;
}

sub _write_fam_file {

## Function : Write pedigree info to file in fam format
## Returns  :
## Arguments: $execution_mode => Either system (direct) or via sbatch
##          : $fam_file_path  => Case file path
##          : $fam_lines_ref  => Fam lines to write
##          : $filehandle     => Filehandle to write to {Optional unless execution_mode=sbatch}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fam_file_path;
    my $fam_lines_ref;
    my $filehandle;

    ## Default(s)
    my $execution_mode;

    my $tmpl = {
        execution_mode => {
            allow       => [qw{ sbatch system }],
            default     => q{sbatch},
            store       => \$execution_mode,
            strict_type => 1,
        },
        fam_lines_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$fam_lines_ref,
            strict_type => 1,
        },
        fam_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$fam_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %execution_mode_map = (
        system => {
            method   => \&_write_fam_file_to_system,
            arg_href => {
                fam_file_path => $fam_file_path,
                fam_lines_ref => $fam_lines_ref,
            },
        },
        sbatch => {
            method   => \&_write_fam_file_to_sbatch,
            arg_href => {
                fam_file_path => $fam_file_path,
                fam_lines_ref => $fam_lines_ref,
                filehandle    => $filehandle,
            },
        },
    );

    $execution_mode_map{$execution_mode}{method}
      ->( { %{ $execution_mode_map{$execution_mode}{arg_href} } } );

    return;
}

sub _write_fam_file_to_sbatch {

## Function : Write pedigree info to file in fam format to sbatch file
## Returns  :
## Arguments: $fam_file_path => Case file path
##          : $fam_lines_ref => Fam lines to write
##          : $filehandle    => Filehandle to write to {Optional unless execution_mode=sbatch}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fam_file_path;
    my $fam_lines_ref;
    my $filehandle;

    my $tmpl = {
        fam_lines_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$fam_lines_ref,
            strict_type => 1,
        },
        fam_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$fam_file_path,
            strict_type => 1,
        },
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Coreutils qw{ gnu_echo };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    say {$filehandle} q{#Generating '.fam' file};

    ## Get parameters
    my @strings = map { $_ . q{\n} } @{$fam_lines_ref};
    gnu_echo(
        {
            enable_interpretation => 1,
            filehandle            => $filehandle,
            no_trailing_newline   => 1,
            outfile_path          => $fam_file_path,
            strings_ref           => \@strings,
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

sub _write_fam_file_to_system {

## Function : Write pedigree info to file in fam format to system
## Returns  :
## Arguments: $fam_file_path => Case file path
##          : $fam_lines_ref => Fam lines to write

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fam_file_path;
    my $fam_lines_ref;

    my $tmpl = {
        fam_lines_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$fam_lines_ref,
            strict_type => 1,
        },
        fam_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$fam_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    # Create anonymous filehandle
    my $filehandle_sys = IO::Handle->new();

    ## Create dir if it does not exists
    make_path( dirname($fam_file_path) );

    open $filehandle_sys, q{>}, $fam_file_path
      or $log->logdie(qq{Can't open $fam_file_path: $ERRNO });

    ## Adds the information from the samples in fam_lines, separated by \n
  LINE:
    foreach my $line ( @{$fam_lines_ref} ) {

        say {$filehandle_sys} $line;
    }
    $log->info( q{Wrote: } . $fam_file_path, $NEWLINE );
    close $filehandle_sys;

    return;
}

1;
