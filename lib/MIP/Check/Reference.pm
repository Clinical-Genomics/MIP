package MIP::Check::Reference;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;
use List::MoreUtils qw { uniq };

##MIPs lib/
use MIP::Recipes::Vt_core qw{ analysis_vt_core };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_references_for_vt check_if_processed_by_vt };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub check_references_for_vt {

##check_references_for_vt

##Function : Check if vt has processed references
##Returns  : @to_process_references
##Arguments: $parameter_href, $active_parameter_href, $vt_references_ref
##         : $parameter_href        => Parameter hash {REF}
##         : $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $vt_references_ref     => The references to check with vt {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $vt_references_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        vt_references_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$vt_references_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Checked references
    my @checked_references;
    ## Store references to process later
    my @to_process_references;

    ## Avoid checking the same reference multiple times
    my %seen;

  PARAMETER_NAME:
    foreach my $parameter_name ( @{$vt_references_ref} ) {

      ASSOCIATED_PROGRAM:
        foreach my $associated_program (
            @{ $parameter_href->{$parameter_name}{associated_program} } )
        {

            ## Alias
            my $active_associated_program =
              $active_parameter_href->{$associated_program};

            next ASSOCIATED_PROGRAM if ( not $active_associated_program );

            ## If SCALAR data type
            if ( $parameter_href->{$parameter_name}{data_type} eq q{SCALAR} ) {

                my $annotation_file = $active_parameter_href->{$parameter_name};

                if ( not exists $seen{$annotation_file} ) {

                    ## Check if vt has processed references using regexp
                    @checked_references = check_if_processed_by_vt(
                        {
                            reference_file_path => $annotation_file,
                            log                 => $log,
                        }
                    );
                    push @to_process_references, @checked_references;
                }
                $seen{$annotation_file} = undef;
            }
            elsif ( $parameter_href->{$parameter_name}{data_type} eq q{ARRAY} )
            {
                ## ARRAY reference

              ANNOTION_FILE:
                foreach my $annotation_file (
                    @{ $active_parameter_href->{$parameter_name} } )
                {

                    if ( not exists $seen{$annotation_file} ) {

                        ## Check if vt has processed references using regexp
                        @checked_references = check_if_processed_by_vt(
                            {
                                reference_file_path => $annotation_file,
                                log                 => $log,
                            }
                        );
                    }
                    push @to_process_references, @checked_references;
                    $seen{$annotation_file} = undef;
                }
            }
            elsif ( $parameter_href->{$parameter_name}{data_type} eq q{HASH} ) {
                ## Hash reference

              ANNOTATION_FILE:
                for my $annotation_file (
                    keys %{ $active_parameter_href->{$parameter_name} } )
                {

                    if ( not exists $seen{$annotation_file} ) {

                        ## Check if vt has processed references using regexp
                        @checked_references = check_if_processed_by_vt(
                            {
                                reference_file_path => $annotation_file,
                                log                 => $log,
                            }
                        );
                    }
                    push @to_process_references, @checked_references;
                    $seen{$annotation_file} = undef;
                }
            }
        }
    }
    return uniq(@to_process_references);
}

sub check_if_processed_by_vt {

##check_if_processed_by_vt

##Function : Check if vt has processed references using regexp
##Returns  : @process_references
##Arguments: $reference_file_path, $log
##         : $reference_file_path => The reference file path
##         : $log                 => Log object
    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_file_path;
    my $log;

    my $tmpl = {
        reference_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$reference_file_path
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %vt_regexp;

    $vt_regexp{vt_decompose}{vcf_key} = q{OLD_MULTIALLELIC};
    $vt_regexp{vt_normalize}{vcf_key} = q{OLD_VARIANT};

    my @to_process_references;
    ## Downloaded and vt later (for downloadable references otherwise
    ## file existens error is thrown downstream)
    if ( not -e $reference_file_path ) {

        ## Do nothing since there is not ref file to check
        return;
    }

  VT_PARAMETER_NAME:
    foreach my $vt_parameter_name ( keys %vt_regexp ) {
        ## MIP flags

        my $regexp =
            q?perl -nae 'if($_=~/ID\=?
          . $vt_regexp{$vt_parameter_name}{vcf_key}
          . q?/) {print $_} if($_=~/#CHROM/) {last}'?;

        ## Detect if vt program has processed reference
        my $ret = `less $reference_file_path | $regexp`;

        ## No trace of vt processing found
        if ( not $ret ) {

            ## Add reference for downstream processing
            push @to_process_references, $reference_file_path;
            $log->warn( q{Cannot detect that }
                  . $vt_parameter_name
                  . q{ has processed reference: }
                  . $reference_file_path
                  . $NEWLINE );
        }
        else {
            ## Found vt processing trace

            $log->info( q{Reference check: }
                  . $reference_file_path
                  . q{ vt: }
                  . $vt_parameter_name
                  . q{ - PASS}
                  . $NEWLINE );
        }
    }
    return uniq(@to_process_references);
}

1;
