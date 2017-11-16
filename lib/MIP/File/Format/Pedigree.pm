package MIP::File::Format::Pedigree;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ make_path };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error};
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Readonly;

## MIPs lib/
use MIP::Gnu::Coreutils qw{ gnu_echo };

BEGIN {
    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ create_fam_file gatk_pedigree_flag reload_previous_pedigree_info };
}

## Constants
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $NEWLINE      => qq{\n};
Readonly my $QUOTE        => q{'};
Readonly my $SPACE        => q{ };
Readonly my $TAB          => qq{\t};
Readonly my $UNDERSCORE   => q{_};

sub create_fam_file {

## Function : Create .fam file to be used in variant calling analyses. Also checks if file already exists when using execution_mode=sbatch.
## Returns  :
## Arguments: $parameter_href        => Hash with paremters from yaml file {REF}
##          : $active_parameter_href => The ac:tive parameters for this analysis hash {REF}
##          : $sample_info_href      => Info on samples and family hash {REF}
##          : $execution_mode        => Either system (direct) or via sbatch
##          : $fam_file_path         => The family file path
##          : $include_header        => Wether to include header ("1") or not ("0")
##          : $FILEHANDLE            => Filehandle to write to {Optional unless execution_mode=sbatch}
##          : $family_id_ref         => The family_id {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $fam_file_path;
    my $FILEHANDLE;

    ## Default(s)
    my $family_id_ref;
    my $execution_mode;
    my $include_header;

    my $tmpl = {
        parameter_href => {
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
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        fam_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$fam_file_path
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
        execution_mode => {
            default     => q{sbatch},
            allow       => [qw{sbatch system}],
            strict_type => 1,
            store       => \$execution_mode
        },
        include_header => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$include_header
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @fam_headers =
      ( q{#family_id}, qw{ sample_id father mother sex phenotype } );

    my @pedigree_lines;
    my $header;

    ## Add @fam_headers
    if ($include_header) {

        push @pedigree_lines, join $TAB, @fam_headers;
    }

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        my $sample_line = ${$family_id_ref};

      HEADER:
        foreach my $header (@fam_headers) {

            if (
                defined $parameter_href->{dynamic_parameter}{$sample_id}
                { q{plink_} . $header } )
            {

                $sample_line .=
                    $TAB
                  . $parameter_href->{dynamic_parameter}{$sample_id}
                  { q{plink_} . $header };
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
                        strings_ref           => \@strings,
                        outfile_path          => $fam_file_path,
                        enable_interpretation => 1,
                        no_trailing_newline   => 1,
                        FILEHANDLE            => $FILEHANDLE,
                    }
                );
                say {$FILEHANDLE} $NEWLINE;
            }
            else {

                $log->fatal(
q{Create fam file[subroutine]:Using 'execution_mode=sbatch' requires a }
                      . q{filehandle to write to. Please supply filehandle to subroutine call}
                      . $NEWLINE );
                exit 1;
            }
        }
    }

    ## Add newly created family file to qc_sample_info
    $sample_info_href->{pedigree_minimal} = $fam_file_path;

    return;
}

sub gatk_pedigree_flag {

## Function : Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
## Returns  : %command
## Arguments: $fam_file_path            => The family file path
##          : $program_name             => The program to use the pedigree file
##          : $pedigree_validation_type => The pedigree validation strictness level

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fam_file_path;
    my $program_name;

    ## Default(s)
    my $pedigree_validation_type;

    my $tmpl = {
        fam_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$fam_file_path
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ SILENT STRICT }],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $parent_counter;
    my $pq_parent_counter = _build_parent_child_counter_regexp(
        {
            family_member => q{parent},
        }
    );

    my $child_counter;
    my $pq_child_counter = _build_parent_child_counter_regexp(
        {
            family_member => q{child},
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

sub reload_previous_pedigree_info {

## Function : Updates sample_info hash with previous run pedigree info
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $sample_info_href      => Info on samples and family hash {REF}
##          : $sample_info_file_path => Previuos sample info file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;
    my $sample_info_file_path;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        sample_info_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_info_file_path
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Yaml qw{ load_yaml };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    if ( -f $active_parameter_href->{sample_info_file} ) {

        $log->info( q{Loaded: } . $active_parameter_href->{sample_info_file},
            $NEWLINE );

        ## Loads a YAML file into an arbitrary hash and returns it.
        # Load parameters from previous run from sample_info_file
        my %previous_sample_info = load_yaml(
            {
                yaml_file => $active_parameter_href->{sample_info_file},
            }
        );

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
## Arguments: $family_member => Parent or child

    my ($arg_href) = @_;

    ## Flatten argument
    my $family_member;

    my $tmpl = {
        family_member => {
            required => 1,
            defined  => 1,
            store    => \$family_member,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    my $regexp = q?perl -ne '?;

    ## Add parent or child, according to which one needs to be counted
    $regexp .= q?my $? . $family_member . $UNDERSCORE . q?counter=0; ?;

    ## Split line around tab if it's not a comment
    $regexp .=
q?while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] eq 0) || ($line[3] eq 0) ) ?;

    ## Increment the counter
    $regexp .= q?{ $?
      . $family_member
      . q?++} } } print $?
      . $family_member
      . q?; last;'?;

    return $regexp;
}

sub _update_sample_info_hash {

## Function : Update sample_info with information from pedigree from previous run. Required e.g. if only updating single sample analysis chains from trio.
## Returns  : %{$previous_sample_info_href}
## Arguments: $previous_sample_info_href => Allowed parameters from pedigre file hash {REF}
##          : $sample_info_href          => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $previous_sample_info_href;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        previous_sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$previous_sample_info_href
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
                my $previous_pedigree_value =
                  delete $sample_href->{$pedigree_key};

                ## Update previous sample info key
                $previous_sample_href->{$pedigree_key} =
                  $previous_pedigree_value;
            }
            else {

                ## New sample_id or key
                $previous_sample_href->{$pedigree_key} =
                  $sample_href->{$pedigree_key};
            }
        }
    }
    return %{$previous_sample_info_href};
}

1;
