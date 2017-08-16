package MIP::File::Format::Pedigree;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case
use Readonly;
use English qw( -no_match_vars );

## MIPs lib/
use MIP::Unix::Standard_streams qw(unix_standard_streams);
use MIP::Unix::Write_to_file qw(unix_write_to_file);
use MIP::Gnu::Coreutils qw (gnu_echo);

BEGIN {
    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(create_fam_file);
}

## Constants
Readonly my $TAB          => qq{\t};
Readonly my $NEWLINE      => qq{\n};
Readonly my $QUOTE        => q{'};
Readonly my $DOUBLE_QUOTE => q{"};

sub create_fam_file {

##create_fam_file

##Function : Create .fam file to be used in variant calling analyses. Also checks if file already exists when using execution_mode=sbatch.
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, sample_info_href, $execution_mode, $fam_file_path, $include_header, $FILEHANDLE, $family_id_ref
##         : $parameter_href        => Hash with paremters from yaml file {REF}
##         : $active_parameter_href => The active parameters for this analysis hash {REF}
##         : $sample_info_href      => Info on samples and family hash {REF}
##         : $execution_mode        => Either system (direct) or via sbatch
##         : $fam_file_path         => The family file path
##         : $include_header        => Wether to include header ("1") or not ("0")
##         : $FILEHANDLE            => Filehandle to write to {Optional unless execution_mode=sbatch}
##         : $family_id_ref         => The family_id {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $execution_mode;
    my $include_header;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $fam_file_path;
    my $FILEHANDLE;

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
            default     => 'sbatch',
            allow       => [ 'sbatch', 'system' ],
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

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger('MIP');

    my @fam_headers =
      ( '#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype' );

    my @pedigree_lines;
    my $header;

    # Add @fam_headers
    if ($include_header) {
        push @pedigree_lines, join $TAB, @fam_headers;
    }

    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {
        my $sample_line = ${$family_id_ref};
        foreach my $header (@fam_headers) {

            if (
                defined(
                    $parameter_href->{dynamic_parameter}{$sample_id}
                      { 'plink_' . $header }
                )
              )
            {
                $sample_line .=
                    $TAB
                  . $parameter_href->{dynamic_parameter}{$sample_id}
                  { 'plink_' . $header };
            }

            elsif (
                defined( $sample_info_href->{sample}{$sample_id}{$header} ) )
            {
                $sample_line .=
                  $TAB . $sample_info_href->{sample}{$sample_id}{$header};
            }
        }
        push @pedigree_lines, $sample_line;
    }

    if ( $execution_mode eq 'system' ) {    #Execute directly
        my $FILEHANDLE_SYS = IO::Handle->new();    #Create anonymous filehandle
        open $FILEHANDLE_SYS, '>', $fam_file_path
          or $log->logdie(qq{Can't open $fam_file_path: $ERRNO $NEWLINE});

     # Adds the information from the samples in pedigree_lines, separatted by \n
        foreach my $line (@pedigree_lines) {
            say {$FILEHANDLE_SYS} $line;
        }
        $log->info( 'Wrote: ' . $fam_file_path, $NEWLINE );
        close $FILEHANDLE_SYS;
    }

    if ( $execution_mode eq 'sbatch' ) {
        if ( !-f $fam_file_path ) {    #Check to see if file already exists
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

1;
