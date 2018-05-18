package MIP::Parse::Parameter;

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
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ parse_infiles parse_prioritize_variant_callers parse_start_with_program };

}

sub parse_infiles {

## Function : Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any files exist.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $indir_path_href       => Indirectories path(s) hash {REF}
##          : $infile_href           => Infiles hash {REF}
##          : $log                   => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $indir_path_href;
    my $infile_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        indir_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$indir_path_href,
            strict_type => 1,
        },
        infile_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_files get_matching_values_key };

    $log->info(q{Reads from platform:});

    ## Collects inputfiles governed by sample_ids
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Return the key if the hash value and query match
        my $infile_directory = get_matching_values_key(
            {
                active_parameter_href => $active_parameter_href,
                parameter_name        => q{infile_dirs},
                query_value           => $sample_id,
            }
        );

        my @infiles = get_files(
            {
                file_directory   => $infile_directory,
                rule_name        => q{*.fastq*},
                rule_skip_subdir => q{original_fastq_files},
            }
        );

        #No "*.fastq*" infiles
        if ( !@infiles ) {

            $log->fatal(
q{Could not find any '.fastq' files in supplied infiles directory }
                  . $infile_directory,
            );
            exit 1;
        }

        ## Check that inFileDirs/infile contains sample_id in filename
      INFILE:
        foreach my $infile (@infiles) {

            unless ( $infile =~ /$sample_id/ ) {

                $log->fatal(
                        q{Could not detect sample_id: }
                      . $sample_id
                      . q{ in supplied infile: }
                      . $infile_directory . q{/}
                      . $infile,
                );
                $log->fatal(
q{Check that: '--sample_ids' and '--inFileDirs' contain the same sample_id and that the filename of the infile contains the sample_id.},
                );
                exit 1;
            }
        }

        $log->info( q{Sample id: } . $sample_id );
        $log->info(qq{\tInputfiles:});

        ## Log each file from platform
      FILE:
        foreach my $file (@infiles) {

            # Indent for visability
            $log->info( qq{\t\t}, $file );
        }

        #Catch inputdir path
        $indir_path_href->{$sample_id} = $infile_directory;

        ## Reload files into hash
        $infile_href->{$sample_id} = [@infiles];
    }
    return;
}

sub parse_prioritize_variant_callers {

## Function : Check that all active variant callers have a prioritization order and that the prioritization elements match a supported variant caller
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

    use MIP::Check::Parameter qw{ check_prioritize_variant_callers };

    my %priority_call_parameter = (
        variant_callers            => q{gatk_combinevariants_prioritize_caller},
        structural_variant_callers => q{sv_svdb_merge_prioritize},
    );

    while ( my ( $variant_caller_type, $prioritize_parameter_name ) =
        each %priority_call_parameter )
    {

        ## Check if we have any active callers
        my $activate_caller_tracker = 0;

      CALLER:
        foreach my $variant_caller (
            @{ $parameter_href->{dynamic_parameter}{$variant_caller_type} } )
        {

            if ( $active_parameter_href->{$variant_caller} ) {

                $activate_caller_tracker++;
            }
        }
        if ($activate_caller_tracker) {

            check_prioritize_variant_callers(
                {
                    active_parameter_href => $active_parameter_href,
                    log                   => $log,
                    parameter_href        => $parameter_href,
                    parameter_name        => $prioritize_parameter_name,
                    variant_callers_ref   => \@{
                        $parameter_href->{dynamic_parameter}
                          {$variant_caller_type}
                    },
                }
            );
        }
        else {
            ## No active callers at all
            return;
        }
    }
    return 1;
}

sub parse_start_with_program {

## Function : Get initiation program, downstream dependencies and update program modes fo start_with_program parameter
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $initiation_file       => Initiation file for pipeline
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $initiation_file;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        initiation_file => {
            defined     => 1,
            required    => 1,
            store       => \$initiation_file,
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

    use MIP::File::Format::Yaml qw{ load_yaml };
    use MIP::Get::Analysis qw{ get_dependency_tree };
    use MIP::Update::Programs qw{  update_program_mode_with_start_with };

    return if ( not defined $active_parameter_href->{start_with_program} );

    my %dependency_tree = load_yaml(
        {
            yaml_file => $initiation_file,
        }
    );

    my @start_with_programs;
    my $is_program_found = 0;
    my $is_chain_found   = 0;

    ## Collects all downstream programs from initation point
    get_dependency_tree(
        {
            dependency_tree_href => \%dependency_tree,
            is_program_found_ref => \$is_program_found,
            is_chain_found_ref   => \$is_chain_found,
            program => $active_parameter_href->{start_with_program},
            start_with_programs_ref => \@start_with_programs,
        }
    );

    ## Update program mode depending on start with flag
    update_program_mode_with_start_with(
        {
            active_parameter_href => $active_parameter_href,
            programs_ref => \@{ $parameter_href->{dynamic_parameter}{program} },
            start_with_programs_ref => \@start_with_programs,
        }
    );
    return 1;
}

1;
