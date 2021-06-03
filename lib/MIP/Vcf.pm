package MIP::Vcf;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $COLON $EQUALS $LOG_NAME $PIPE $SEMICOLON $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_sample_ids_from_vcf };
}

sub get_sample_ids_from_vcf {

## Function : Get sample ids from a vcf file
## Returns  : @sample_ids
## Arguments: $vcf_file_path => Unannotated case vcf file from dna pipeline

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vcf_file_path;

    my $tmpl = {
        vcf_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$vcf_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Program::Bcftools qw{ bcftools_view };
    use MIP::Program::Gnu::Bash qw{ gnu_export gnu_unset };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Export MIP_BIND to bind reference path to htslib sif in proxy bin
    my @commands = gnu_export(
        {
            bash_variable => q{MIP_BIND} . $EQUALS . $vcf_file_path . $COLON . $vcf_file_path,
        }
    );
    push @commands, $SEMICOLON;

    ## Get sample_ids from vcf
    push @commands,
      bcftools_view(
        {
            header_only => 1,
            infile_path => $vcf_file_path,
        }
      );

    push @commands, $PIPE;

    push @commands,
      perl_nae_oneliners(
        {
            oneliner_name => q{get_vcf_sample_ids},
        }
      );
    push @commands, $SEMICOLON;

    ## Unset MIP_BIND after system parsing
    push @commands, gnu_unset( { bash_variable => q{MIP_BIND}, } );

    my %process_return = child_process(
        {
            commands_ref => \@commands,
            process_type => q{open3},
        }
    );

    if ( @{ $process_return{stderrs_ref} }
        or not @{ $process_return{stdouts_ref} } )
    {

        $log->fatal(qq{Could not retrieve sample id from vcf: $vcf_file_path});

      ERROR_LINE:
        foreach my $error_line ( @{ $process_return{stderrs_ref} } ) {
            $log->fatal(qq{ERROR: $error_line});
        }
        exit 1;
    }

    my @sample_ids = split $SPACE, $process_return{stdouts_ref}->[0];

    return @sample_ids;
}

1;
