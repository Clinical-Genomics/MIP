package MIP::File::Format::Pedigree;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ make_path };
use List::Util qw{ none };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error};
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use List::MoreUtils qw { any };
use Readonly;

## MIPs lib/
use MIP::Gnu::Coreutils qw{ gnu_echo };
use MIP::Constants qw{ $NEWLINE $SPACE $TAB $UNDERSCORE };

BEGIN {
    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.11;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ create_fam_file detect_founders detect_sample_id_gender detect_trio gatk_pedigree_flag has_trio is_sample_proband_in_trio parse_yaml_pedigree_file reload_previous_pedigree_info };
}

## Constants
Readonly my $TRIO_MEMBERS_COUNT => 3;

sub is_sample_proband_in_trio {

## Function : Check if sample id has an affected or unknown phenotype and is child in trio
## Returns  : 0 | 1
## Arguments: $sample_id        => Sample id
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
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

    use MIP::Get::Parameter qw{ get_pedigree_sample_id_attributes };
    use MIP::Sample_info qw{ get_family_member_id };

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
    return 0 if ( $phenotype eq q{unaffected} );

    ## Get family hash
    my %family_member_id =
      get_family_member_id( { sample_info_href => $sample_info_href } );

    ## Check if the sample is an affected child
    return 0 if ( none { $_ eq $sample_id } @{ $family_member_id{children} } );

    return 1;
}

sub create_fam_file {

## Function : Create .fam file to be used in variant calling analyses. Also checks if file already exists when using execution_mode=sbatch.
## Returns  :
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $case_id               => Case_id
##          : $execution_mode        => Either system (direct) or via sbatch
##          : $fam_file_path         => Case file path
##          : $FILEHANDLE            => Filehandle to write to {Optional unless execution_mode=sbatch}
##          : $include_header        => Include header ("1") or not ("0")
##          : $log                   => Log object
##          : $parameter_href        => Hash with paremters from yaml file {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $fam_file_path;
    my $FILEHANDLE;
    my $log;
    my $parameter_href;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $execution_mode;
    my $include_header;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        execution_mode => {
            allow       => [qw{sbatch system}],
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
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        include_header => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$include_header,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
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

    my @fam_headers = ( q{#family_id}, qw{ sample_id father mother sex phenotype } );

    my @pedigree_lines;
    my $header;

    ## Add @fam_headers
    if ($include_header) {

        push @pedigree_lines, join $TAB, @fam_headers;
    }

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        my $sample_line = $case_id;

      HEADER:
        foreach my $header (@fam_headers) {

            if ( defined $parameter_href->{cache}{$sample_id}{ q{plink_} . $header } ) {

                $sample_line .=
                  $TAB . $parameter_href->{cache}{$sample_id}{ q{plink_} . $header };
            }
            elsif ( defined $sample_info_href->{sample}{$sample_id}{$header} ) {

                $sample_line .=
                  $TAB . $sample_info_href->{sample}{$sample_id}{$header};
            }
        }
        push @pedigree_lines, $sample_line;
    }

    ## Execute directly
    if ( $execution_mode eq q{system} ) {

        # Create anonymous filehandle
        my $FILEHANDLE_SYS = IO::Handle->new();

        ## Create dir if it does not exists
        make_path( dirname($fam_file_path) );

        open $FILEHANDLE_SYS, q{>}, $fam_file_path
          or $log->logdie(qq{Can't open $fam_file_path: $ERRNO });

        ## Adds the information from the samples in pedigree_lines, separated by \n
      LINE:
        foreach my $line (@pedigree_lines) {

            say {$FILEHANDLE_SYS} $line;
        }
        $log->info( q{Wrote: } . $fam_file_path, $NEWLINE );
        close $FILEHANDLE_SYS;
    }

    if ( $execution_mode eq q{sbatch} ) {

        ## Check to see if file already exists
        if ( not -f $fam_file_path ) {

            if ($FILEHANDLE) {

                say {$FILEHANDLE} q{#Generating '.fam' file};

                ## Get parameters
                my @strings = map { $_ . q{\n} } @pedigree_lines;
                gnu_echo(
                    {
                        enable_interpretation => 1,
                        FILEHANDLE            => $FILEHANDLE,
                        no_trailing_newline   => 1,
                        outfile_path          => $fam_file_path,
                        strings_ref           => \@strings,
                    }
                );
                say {$FILEHANDLE} $NEWLINE;
            }
            else {

                $log->fatal(
q{Create fam file[subroutine]:Using 'execution_mode=sbatch' requires a }
                      . q{filehandle to write to. Please supply filehandle to subroutine call}
                );
                exit 1;
            }
        }
    }

    ## Add newly created case file to qc_sample_info
    $sample_info_href->{pedigree_minimal} = $fam_file_path;

    return;
}

sub detect_founders {

## Function  : Detect number of founders (i.e. parents ) based on pedigree file
## Returns   : "scalar @founders"
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

    my @founders;

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Alias
        my $father_id =
          $sample_info_href->{sample}{$sample_id}{father};
        my $mother_id =
          $sample_info_href->{sample}{$sample_id}{mother};

        my %parent = (
            father => $father_id,
            mother => $mother_id,
        );

      PARENT:
        foreach my $parent_id ( values %parent ) {

            ## Child for true parent_id
            next PARENT if ( not defined $parent_id and not $parent_id );

            ## If parent is present in current analysis
            if (
                any { $_ eq $parent_id }
                @{ $active_parameter_href->{sample_ids} }
              )
            {

                push @founders, $parent_id;
            }
        }
    }
    return scalar @founders;
}

sub detect_sample_id_gender {

## Function : Detect gender of the current analysis
## Returns  : "$found_male $found_female $found_other $found_other_count"
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

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

    ## Initialize
    my $found_male        = 0;
    my $found_female      = 0;
    my $found_other       = 0;
    my $found_other_count = 0;

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## If male
        if ( $sample_info_href->{sample}{$sample_id}{sex} =~ / 1 | ^male/sxm ) {

            $found_male = 1;
        }
        elsif ( $sample_info_href->{sample}{$sample_id}{sex} =~ / 2 | female /sxm ) {
            ## If female

            $found_female = 1;
        }
        else {
            ## Must be other

            ## Include since it might be male to enable analysis of Y.
            $found_male = 1;

            # "Other" metrics
            $found_other = 1;
            $found_other_count++;
        }
    }
    return $found_male, $found_female, $found_other, $found_other_count;
}

sub detect_trio {

## Function  : Detect case constellation based on pedigree file
## Returns   : undef | 1
## Arguments : $active_parameter_href => Active parameters for this analysis hash {REF}
##           : $log                   => Log
##           : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $sample_info_href;

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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( scalar @{ $active_parameter_href->{sample_ids} } == 1 ) {

        $log->info( q{Found single sample: } . $active_parameter_href->{sample_ids}[0] );
        return;
    }
    elsif ( scalar @{ $active_parameter_href->{sample_ids} } == $TRIO_MEMBERS_COUNT ) {

        my $is_trio;

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Alias
            my $father_id =
              $sample_info_href->{sample}{$sample_id}{father};
            my $mother_id =
              $sample_info_href->{sample}{$sample_id}{mother};
            $is_trio = _parse_trio_members(
                {
                    father_id      => $father_id,
                    log            => $log,
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

    use MIP::Get::Parameter qw{ get_pedigree_sample_id_attributes };

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

sub parse_yaml_pedigree_file {

## Function : Parse case info in YAML pedigree file. Check pedigree data for allowed entries and correct format. Add data to sample_info and active_parameter depending on user info.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_path             => Pedigree file path
##          : $log                   => Log object
##          : $parameter_href        => Parameter hash {REF}
##          : $pedigree_href         => Pedigree hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_path;
    my $log;
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
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
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

    use MIP::Check::Pedigree
      qw{ check_founder_id check_pedigree_mandatory_key check_pedigree_sample_allowed_values check_pedigree_vs_user_input_sample_ids };
    use MIP::Get::Parameter qw{ get_capture_kit get_user_supplied_info };
    use MIP::Set::Pedigree
      qw{ set_active_parameter_pedigree_keys set_pedigree_capture_kit_info set_pedigree_case_info set_pedigree_phenotype_info set_pedigree_sample_info set_pedigree_sex_info };

    ## Use to collect which sample_ids have used a certain capture_kit
    my $case_id = $pedigree_href->{case};

    ## User supplied sample_ids via cmd or config
    my @user_input_sample_ids;

    ## Check pedigree mandatory keys
    check_pedigree_mandatory_key(
        {
            file_path     => $file_path,
            log           => $log,
            pedigree_href => $pedigree_href,
        }
    );

    ## Check that supplied cmd and YAML pedigree case_id match
    if ( $pedigree_href->{case} ne $active_parameter_href->{case_id} ) {

        $log->fatal( q{Pedigree file: }
              . $file_path
              . q{ for  pedigree case_id: '}
              . $pedigree_href->{case}
              . q{' and supplied case: '}
              . $active_parameter_href->{case_id}
              . q{' does not match} );
        exit 1;
    }

    ### Check sample keys values
    check_pedigree_sample_allowed_values(
        {
            file_path     => $file_path,
            log           => $log,
            pedigree_href => $pedigree_href,
        }
    );

    ## Get potential input from user
    my %user_supply_switch = get_user_supplied_info(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    if ( $user_supply_switch{sample_ids}
        and exists $active_parameter_href->{sample_ids} )
    {

        ## Set cmd or config supplied sample_ids
        @user_input_sample_ids = @{ $active_parameter_href->{sample_ids} };
    }

    set_pedigree_case_info(
        {
            pedigree_href    => $pedigree_href,
            sample_info_href => $sample_info_href,
        }
    );

    my @pedigree_sample_ids = set_pedigree_sample_info(
        {
            active_parameter_href     => $active_parameter_href,
            pedigree_href             => $pedigree_href,
            sample_info_href          => $sample_info_href,
            user_supply_switch_href   => \%user_supply_switch,
            user_input_sample_ids_ref => \@user_input_sample_ids,
        }
    );

    ## Add sex to dynamic parameters
    set_pedigree_sex_info(
        {
            pedigree_href  => $pedigree_href,
            parameter_href => $parameter_href,
        }
    );

    ## Add phenotype to dynamic parameters
    set_pedigree_phenotype_info(
        {
            pedigree_href  => $pedigree_href,
            parameter_href => $parameter_href,
        }
    );

    set_active_parameter_pedigree_keys(
        {
            active_parameter_href   => $active_parameter_href,
            pedigree_href           => $pedigree_href,
            sample_info_href        => $sample_info_href,
            user_supply_switch_href => \%user_supply_switch,
        }
    );

    set_pedigree_capture_kit_info(
        {
            active_parameter_href   => $active_parameter_href,
            parameter_href          => $parameter_href,
            pedigree_href           => $pedigree_href,
            sample_info_href        => $sample_info_href,
            user_supply_switch_href => \%user_supply_switch,
        }
    );

    ## Check that founder_ids are included in the pedigree info and the analysis run
    check_founder_id(
        {
            active_sample_ids_ref => \@{ $active_parameter_href->{sample_ids} },
            log                   => $log,
            pedigree_href         => $pedigree_href,
        }
    );

    if ( not $user_supply_switch{sample_ids} ) {

        ## Lexiographical sort to determine the correct order of ids indata
        @{ $active_parameter_href->{sample_ids} } =
          sort @{ $active_parameter_href->{sample_ids} };
    }
    else {

        ## Check that CLI supplied sample_id exists in pedigree
        check_pedigree_vs_user_input_sample_ids(
            {
                file_path                 => $file_path,
                log                       => $log,
                pedigree_sample_ids_ref   => \@pedigree_sample_ids,
                user_input_sample_ids_ref => \@user_input_sample_ids,
            }
        );
    }
    return 1;
}

sub reload_previous_pedigree_info {

## Function : Updates sample_info hash with previous run pedigree info
## Returns  :
## Arguments: $log                   => Log object to write to
##          : $sample_info_file_path => Previuos sample info file
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $sample_info_file_path;
    my $sample_info_href;

    my $tmpl = {
        log                   => { required => 1, store => \$log, },
        sample_info_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$sample_info_file_path,
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

    use MIP::File::Format::Yaml qw{ load_yaml };

    if ( -f $sample_info_file_path ) {

        ## Loads a YAML file into an arbitrary hash and returns it.
        # Load parameters from previous run from sample_info_file
        my %previous_sample_info = load_yaml(
            {
                yaml_file => $sample_info_file_path,
            }
        );

        $log->info( q{Loaded: } . $sample_info_file_path, $NEWLINE );

        ## Update sample_info with pedigree information from previous run
        ## Should be only pedigree keys in %allowed_entries
        %{$sample_info_href} = _update_sample_info_hash(
            {
                sample_info_href          => $sample_info_href,
                previous_sample_info_href => \%previous_sample_info,
            }
        );
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

sub _update_sample_info_hash {

## Function : Update sample_info with information from pedigree from previous run. Required e.g. if only updating single sample analysis chains from trio.
## Returns  : %{$previous_sample_info_href}
## Arguments: $previous_sample_info_href => Allowed parameters from pedigre file hash {REF}
##          : $sample_info_href          => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $previous_sample_info_href;
    my $sample_info_href;

    my $tmpl = {
        previous_sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$previous_sample_info_href,
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

  SAMPLE_ID:
    foreach my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

        ## Alias
        my $sample_href = \%{ $sample_info_href->{sample}{$sample_id} };
        my $previous_sample_href =
          \%{ $previous_sample_info_href->{sample}{$sample_id} };

      PEDIGREE_KEY:
        foreach my $pedigree_key ( keys %{$sample_href} ) {

            ## Previous run information, which should be updated using pedigree from current analysis
            if ( exists $previous_sample_href->{$pedigree_key} ) {

                ## Required to update keys downstream
                my $previous_pedigree_value = delete $sample_href->{$pedigree_key};

                ## Update previous sample info key
                $previous_sample_href->{$pedigree_key} = $previous_pedigree_value;
            }
            else {

                ## New sample_id or key
                $previous_sample_href->{$pedigree_key} = $sample_href->{$pedigree_key};
            }
        }
    }
    return %{$previous_sample_info_href};
}

sub _parse_trio_members {

## Function  : Parse trio constellation
## Returns   : %trio
## Arguments : $father_id      => Potential father
##           : $log            => Log
##           : $mother_id      => Potential mother
##           : $sample_id      => Sample under investigation
##           : $sample_ids_ref => Sample_ids in current analysis {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $father_id;
    my $log;
    my $mother_id;
    my $sample_id;
    my $sample_ids_ref;

    my $tmpl = {
        father_id => {
            required => 1,
            store    => \$father_id,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
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
        if (
            any { $_ eq $parent_id }
            @{$sample_ids_ref}
          )
        {

            ## Set as parents
            $trio{$parent_role} = $parent_id;
        }
    }
    if ( scalar( keys %trio ) == $TRIO_MEMBERS_COUNT ) {

        $log->info(
                q{Found trio: Child = "}
              . $trio{child}
              . q{", Father = "}
              . $trio{father}
              . q{", Mother = "}
              . $trio{mother} . q{"},
        );
        return 1;
    }
    return;
}

1;
