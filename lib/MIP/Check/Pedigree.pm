package MIP::Check::Pedigree;

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
use List::MoreUtils qw { any };
use Readonly;

## MIP lib
use MIP::Constants qw{ $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_founder_id
      check_pedigree_mandatory_key
      check_pedigree_sample_allowed_values
      check_pedigree_vs_user_input_sample_ids
    };
}

sub check_founder_id {

## Function : Check that founder_ids are included in the pedigree info
## Returns  :
## Arguments: $active_sample_ids_ref => Array of pedigree samples {REF}
##          : $log                   => Log object
##          : $pedigree_href         => Pedigree info {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_sample_ids_ref;
    my $log;
    my $pedigree_href;

    my $tmpl = {
        active_sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$active_sample_ids_ref,
            strict_type => 1,
        },
        log           => { defined => 1, required => 1, store => \$log, },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

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
##          : $file_path             => Pedigree file path
##          : $pedigree_href         => YAML pedigree info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_path;
    my $pedigree_href;

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

    ### Family
    ## Check mandatory case keys
  MANDATORY_KEY:
    foreach my $key (@mandatory_case_keys) {

        next MANDATORY_KEY if ( exists $pedigree_href->{$key} );

        $log->fatal( q{Pedigree file: }
              . $file_path
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
                  . $file_path
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
## Arguments: $file_path     => Pedigree file path
##          : $log           => Log object
##          : $pedigree_href => YAML pedigree info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $log;
    my $pedigree_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        log           => { store => \$log },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %allowed_values = (
        analysis_type => [qw{ dragen_rd_dna vrn wes wgs wts }],
        phenotype     => [qw{ affected unaffected unknown }],
        sex           => [qw{ female male other unknown }],
    );

    ## Check input to sample_info hash for at sample level
  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

      SAMPLE_KEY:
        foreach my $key ( keys %{$pedigree_sample_href} ) {

            ## No defined allowed values
            next SAMPLE_KEY if ( not exists $allowed_values{$key} );

            ## If element is not part of array
            next SAMPLE_KEY
              if (
                any { $_ eq $pedigree_sample_href->{$key} }
                @{ $allowed_values{$key} }
              );

            $log->fatal( q{Pedigree file: } . $file_path );
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
## Arguments: $file_path                 => Pedigree file path
##          : $log                       => Log object
##          : $pedigree_sample_ids_ref   => Pedigree sample ids
##          : $user_input_sample_ids_ref => user input sample ids

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $log;
    my $pedigree_sample_ids_ref;
    my $user_input_sample_ids_ref;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        log                     => { required => 1, store => \$log },
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

  SAMPLE_ID:
    foreach my $sample_id ( @{$user_input_sample_ids_ref} ) {

        ## If element is not part of array
        if ( not any { $_ eq $sample_id } @{$pedigree_sample_ids_ref} ) {

            $log->fatal( q{File: }
                  . $file_path
                  . q{ provided sample_id: }
                  . $sample_id
                  . q{ is not present in pedigree file} );
            exit 1;
        }
    }
    return 1;
}

1;
