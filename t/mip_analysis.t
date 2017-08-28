#!/usr/bin/env perl

###Will test perl modules and some selected funtions as well as vcf keys both in header and body. Adjusts dynamically according to supplied config file.
###Copyright 2016 Henrik Stranneheim

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie qw(open close :all);
use 5.018;    #Require at least perl 5.18
use utf8;     #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);

use Cwd qw(abs_path);
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile devnull);
use FindBin qw($Bin);    #Find directory of script
use Getopt::Long;
use Params::Check qw[check allow last_error];
use Test::More;

## Third party module(s)
use List::Util qw(any);

##MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use File::Format::Yaml qw(load_yaml);
use MIP_log::Log4perl qw(initiate_logger);
use Check::Check_modules qw(check_modules);
use Script::Utils qw(help);

our $USAGE = build_usage( {} );

BEGIN {

    my @modules = (
        'Modern::Perl',              #MIP
        'autodie',                   #MIP
        'IPC::System::Simple',       #MIP
        'Path::Iterator::Rule',      #MIP
        'YAML',                      #MIP
        'File::Format::Yaml',        #MIP
        'Log::Log4perl',             #MIP
        'MIP_log::Log4perl',         #MIP
        'List::Util',                #MIP
        'Set::IntervalTree',         # MIP/vcfParser.pl
        'Net::SSLeay',               # VEP
        'LWP::Simple',               # VEP
        'LWP::Protocol::https',      # VEP
        'PerlIO::gzip',              #VEP
        'IO::Uncompress::Gunzip',    #VEP
        'HTML::Lint',                #VEP
        'Archive::Zip',              # VEP
        'Archive::Extract',          #VEP
        'DBI',                       # VEP
        'JSON',                      # VEP
        'DBD::mysql',                # VEP
        'CGI',                       # VEP
        'Sereal::Encoder',           # VEP
        'Sereal::Decoder',           # VEP
        'Bio::Root::Version',        #VEP
        'Module::Build',             #VEP
        'File::Copy::Recursive',     #VEP
    );

    ## Evaluate that all modules required are installed
    Check::Check_modules::check_modules(
        {
            modules_ref  => \@modules,
            program_name => $PROGRAM_NAME,
        }
    );
}

my ( $infile, $config_file );
my ( %parameter, %active_parameter, %pedigree, %vcfparser_data,
    %consequence_severity );

my $VERSION = '2.0.0';

if ( scalar(@ARGV) == 0 ) {

    print STDOUT $USAGE, "\n";
    exit;
}

############
####MAIN####
############

$infile      = $ARGV[0];
$config_file = $ARGV[1];

###User Options
GetOptions(
    'h|help' => sub { print STDOUT $USAGE, "\n"; exit; },    #Display help text
    'v|version' => sub {
        print STDOUT "\n" . basename($PROGRAM_NAME) . q{ } . $VERSION, "\n\n";
        exit;
    },    #Display version number
  )
  or Script::Utils::help(
    {
        USAGE     => $USAGE,
        exit_code => 1,
    }
  );

unless ( defined($infile) ) {

    print STDERR 'Please supply an infile', "\n";
    exit;
}

unless ( defined($config_file) ) {

    print STDERR 'Please supply a config file', "\n";
    exit;
}

## Test perl modules and functions
test_modules();

if ( defined($config_file) ) {    #Input from cmd

    ## Loads a YAML file into an arbitrary hash and returns it.
    %active_parameter = load_yaml( { yaml_file => $config_file, } );
}

if ( exists( $active_parameter{pedigree_file} ) ) {

    ## Loads a YAML file into an arbitrary hash and returns it.
    %pedigree = load_yaml( { yaml_file => $active_parameter{pedigree_file}, } );
    ### Sample level info
    foreach my $pedigree_sample_href ( @{ $pedigree{samples} } ) {

        ## Sample_id
        my $sample_id = $pedigree_sample_href->{sample_id};    #Alias

        ## Phenotype
        push
          @{ $parameter{dynamic_parameter}{ $pedigree_sample_href->{phenotype} }
          }, $sample_id;

        ## Sex
        push @{ $parameter{dynamic_parameter}{ $pedigree_sample_href->{sex} } },
          $sample_id;
    }
}

if ( $infile =~ /.selected.vcf/xm ) {

    if (   ( defined( $active_parameter{vcfparser_select_file} ) )
        && ( $active_parameter{vcfparser_select_file} ) )
    {

        ## Reads a file containg features to be annotated using range queries
        read_range_file(
            {
                vcfparser_data_href => \%vcfparser_data,
                range_coulumns_ref  => \@{
                    $active_parameter{vcfparser_select_feature_annotation_columns}
                },
                infile_path =>
                  catfile( $active_parameter{vcfparser_select_file} ),
                range_file_key => 'select_file',
            }
        );
    }
}
else {    #Range file

    if (   ( defined( $active_parameter{vcfparser_range_feature_file} ) )
        && ( $active_parameter{vcfparser_range_feature_file} ) )
    {

        ## Reads a file containg features to be annotated using range queries
        read_range_file(
            {
                vcfparser_data_href => \%vcfparser_data,
                range_coulumns_ref  => \@{
                    $active_parameter{vcfparser_range_feature_annotation_columns}
                },
                infile_path =>
                  catfile( $active_parameter{vcfparser_range_feature_file} ),
                range_file_key => 'range_file',
            }
        );
    }
}

## Reads infile in vcf format and parses annotations
read_infile_vcf(
    {
        parameter_href            => \%parameter,
        active_parameter_href     => \%active_parameter,
        vcfparser_data_href       => \%vcfparser_data,
        consequence_severity_href => \%consequence_severity,
        infile                    => $infile,
    }
);

# Reached the end safely
Test::More::done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $program_name
##         : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    return <<"END_USAGE";
 $program_name infile.vcf [VCF] config_file [YAML] [options]
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}

sub test_modules {

##test_modules

##Function : Test perl modules and functions
##Returns  : ""
##Arguments:
##         :

    print STDOUT "\nTesting perl modules and selected functions\n\n";

    use FindBin qw($Bin);    #Find directory of script

    ok( defined($Bin), 'FindBin: Locate directory of script' );

    use File::Basename qw(dirname);    #Strip the last part of directory

    ok( dirname($Bin),
        'File::Basename qw(dirname): Strip the last part of directory' );

    use File::Spec::Functions qw(catdir);

    ok( catdir( dirname($Bin), 't' ),
        'File::Spec::Functions qw(catdir): Concatenate directories' );

    use YAML;

    my $yaml_file =
      catdir( dirname($Bin), qw(definitions define_parameters.yaml) );
    ok( -f $yaml_file, 'YAML: File= ' . $yaml_file . 'in MIP directory' );

    my $yaml = YAML::LoadFile($yaml_file);    #Create an object
    ok( defined $yaml, 'YAML: Load File' );   #Check that we got something
    ok( Dump($yaml),   'YAML: Dump file' );

    use Log::Log4perl;
    ## Creates log
    my $log_file = catdir( dirname($Bin), qw(templates mip_config.yaml) );
    ok( -f $log_file,
        'Log::Log4perl: File= ' . $log_file . 'in MIP directory' );

    ## Creates log object
    my $log = initiate_logger(
        {
            categories_ref => [qw(TRACE ScreenApp)],
            file_path_ref  => \$log_file,
            log_name       => 'Test',
        }
    );

    ok( $log->info(1),  'Log::Log4perl: info' );
    ok( $log->warn(1),  'Log::Log4perl: warn' );
    ok( $log->error(1), 'Log::Log4perl: error' );
    ok( $log->fatal(1), 'Log::Log4perl: fatal' );

    use Getopt::Long;
    push @ARGV, qw(-verbose 2);
    my $verbose = 1;
    ok(
        GetOptions( 'verbose:n' => \$verbose ),
        'Getopt::Long: Get options call'
    );
    ok( $verbose == 2, 'Getopt::Long: Get options modified' );

    return;
}

sub read_infile_vcf {

##read_infile_vcf

##Function : Reads infile in vcf format and adds and parses annotations
##Returns  : ""
##Arguments: $parameter_href, $active_parameter_href, $vcfparser_data_href, $consequence_severity_href, $infile
##         : $parameter_href            => The parameter hash {REF}
##         : $active_parameter_href     => The active parameters for this analysis hash {REF}
##         : vcfparser_data_href        => The keys from vcfParser i.e. range file and select file
##         : $consequence_severity_href => Consequence severity for SO-terms {REF}
##         : $infile                    => Infile

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $vcfparser_data_href;
    my $consequence_severity_href;
    my $infile;

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
        vcfparser_data_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$vcfparser_data_href
        },
        consequence_severity_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$consequence_severity_href
        },
        infile =>
          { required => 1, defined => 1, strict_type => 1, store => \$infile },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger('Test');

    my @vep_format_fields;
    my %vep_format_field_column;

    my %meta_data;
    my %vcf_header;
    my @vcf_format_columns;    #Catch #vcf header #CHROM line
    my %vcf_info_key;
    my %vcf_info_csq_key;

    $log->info( 'Testing vcf header in file: ' . $infile, "\n\n" );

    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle
    open( $FILEHANDLE, '<', $infile )
      or $log->logdie( 'Cannot open ' . $infile . ':' . $!, "\n" );

    while (<$FILEHANDLE>) {

        chomp $_;

        # Exit after parsing X lines
        if ( $. > 10000 ) {
            last;
        }

        # Avoid blank lines
        if (m/^\s+$/) {
            next;
        }

        # MetaData
        if ( $_ =~ /^##(\S+)=/ ) {

            parse_meta_data( \%meta_data, $_ );

            if ( $_ =~ /INFO\=\<ID\=([^,]+)/ ) {

                $vcf_header{INFO}{$1} = $1;
            }

            # Find VEP INFO Field
            if ( $_ =~ /INFO\=\<ID\=CSQ/ ) {

                ok(
                    $active_parameter_href->{pvarianteffectpredictor},
                    'VEP: CSQ in header and VEP should have been executed'
                );

                # Locate Format within VEP INFO meta line
                if ( $_ =~ /Format:\s(\S+)"\>/ ) {

                    @vep_format_fields = split /\|/, $1;

                    while ( my ( $field_index, $field ) =
                        each(@vep_format_fields) )
                    {

                        # Save the order of VEP features
                        $vep_format_field_column{$field} = $field_index;
                    }
                }
                next;
            }
            next;
        }
        if ( $_ =~ /^#CHROM/ ) {

            # Split vcf format line
            @vcf_format_columns = split /\t/, $_;

            ### Check Header now that we read all

            ## VT
            if ( $active_parameter_href->{pvt} > 0 ) {

                if ( $active_parameter_href->{vt_decompose} > 0 ) {

                    ok( defined( $vcf_header{INFO}{OLD_MULTIALLELIC} ),
                        'VTDecompose key: OLD_MULTIALLELIC' );
                }
                if ( $active_parameter_href->{vt_normalize} > 0 ) {

                    ok( defined( $vcf_header{INFO}{OLD_VARIANT} ),
                        'VTNormalize key: OLD_VARIANT' );
                }
            }

            ## VCFParser
            if ( $active_parameter_href->{pvcfparser} > 0 ) {

                for my $key ( keys %{$vcfparser_data_href} ) {

                    ok( defined( $vcf_header{INFO}{$key} ),
                        'VCFParser key: ' . $key );
                }

                ## Keys from vcfParser that are dynamically created from parsing the data
                my @vcfparser_dynamic_keys = ( 'most_severe_consequence', );

                foreach my $key (@vcfparser_dynamic_keys) {

                    ok(
                        defined( $vcf_header{INFO}{$key} ),
                        'vcfparser dynamic keys: ' . $key
                    );
                }
            }

            ## Snpeff
            if ( $active_parameter_href->{psnpeff} > 0 ) {

                my @splitted_keys;
                for my $annotation_file (
                    keys %{ $active_parameter_href->{snpsift_annotation_files} }
                  )
                {

                    if (
                        $active_parameter_href->{snpsift_annotation_outinfo_key}
                        {$annotation_file} )
                    {

                        @splitted_keys = split ',', $active_parameter_href
                          ->{snpsift_annotation_outinfo_key}{$annotation_file};

                        my @original_vcf_keys = split ',',
                          $active_parameter_href->{snpsift_annotation_files}
                          {$annotation_file};

                        ## Modify list elements in place to produce -names flag from SnpEff
                        foreach my $elements (@splitted_keys) {

                            $elements .= shift(@original_vcf_keys);
                        }
                    }
                    else {

                        @splitted_keys = split ',',
                          $active_parameter_href->{snpsift_annotation_files}
                          {$annotation_file};
                    }
                    foreach my $key (@splitted_keys) {

                        ok(
                            defined( $vcf_header{INFO}{$key} ),
                            'Snpsift annotation key: ' . $key
                        );
                    }
                }
                for my $key (
                    @{ $active_parameter_href->{snpsift_dbnsfp_annotations} } )
                {

                    $key =~ s/\+/_/g
                      ; #Special case due to the fact that snpEff v4.2 transforms + to _ for some reason
                    $key = 'dbNSFP_' . $key;
                    ok(
                        defined( $vcf_header{INFO}{$key} ),
                        'Snpsift dbNSFP_key: ' . $key
                    );
                }
            }
            ## GENMOD
            ## 1000G key from genmodFilter
            if (   ( $active_parameter_href->{pvt} > 0 )
                && ( $active_parameter_href->{vt_genmod_filter} > 0 ) )
            {

                if (
                    defined( $active_parameter_href->{vt_genmod_filter_1000G} )
                  )
                {

                    ok( defined( $vcf_header{INFO}{'1000GAF'} ),
                        'Genmod filter: 1000GAF key' );
                }
            }
            if ( $active_parameter_href->{prankvariant} > 0 ) {

                ## Keys from genmod
                my @genmod_keys;

                if (
                    (
                        defined(
                            $parameter_href->{dynamic_parameter}{unaffected}
                        )
                    )
                    && ( @{ $parameter_href->{dynamic_parameter}{unaffected} }
                        eq @{ $active_parameter_href->{sample_ids} } )
                  )
                {    #Only unaffected - do nothing
                }
                else {

                    @genmod_keys = (
                        'Compounds',  'RankScore',
                        'ModelScore', 'GeneticModels',
                    );
                }

                foreach my $key (@genmod_keys) {

                    ok( defined( $vcf_header{INFO}{$key} ), 'Genmod: ' . $key );
                }

                ## Spidex key from genmodAnnotate
                if ( $active_parameter_href->{spidex_file} ) {

                    ok(
                        defined( $vcf_header{INFO}{SPIDEX} ),
                        'Genmod annotate: SPIDEX key'
                    );
                }
                ## CADD key from genmodAnnotate
                if (   ( $active_parameter_href->{genmod_annotate_cadd_files} )
                    || ( $active_parameter_href->{genmod_annotate_cadd_files} )
                  )
                {

                    ok(
                        defined( $vcf_header{INFO}{CADD} ),
                        'Genmod annotate: CADD key'
                    );
                }
            }
            next;
        }
        if ( $_ =~ /^(\S+)/ ) {

            my %record;

            my @line_elements = split "\t", $_;    #Loads vcf elements

            ##Add line elements to record hash
          LINE:
            while ( my ( $element_index, $element ) = each(@line_elements) ) {

                $record{ $vcf_format_columns[$element_index] } =
                  $element;    #Link vcf format headers to the line elements
            }

            my @info_elements = split /;/,
              $record{INFO};    #Split INFO field to key=value items

            ##Add INFO to record hash as separate key
          INFO:
            for my $element (@info_elements) {

                my @key_value_pairs = split "=",
                  $element;     #key index = 0 and value index = 1

                $vcf_info_key{ $key_value_pairs[0] }++;    #Increment
                $record{INFO_key_value}{ $key_value_pairs[0] } =
                  $key_value_pairs[1];
            }
            if ( exists( $record{INFO_key_value}{CSQ} ) ) {

                my @transcripts = split /,/,
                  $record{INFO_key_value}{CSQ};    #Split into transcripts

              CSQ_TRANSCRIPT:
                foreach my $transcript (@transcripts) {

                    my @transcript_effects = split /\|/, $transcript;

                  CSQ_TRANSCRIPT_EFFECTS:
                    foreach my $effect ( keys %vep_format_field_column ) {

                        if (
                            (
                                defined(
                                    $transcript_effects
                                      [ $vep_format_field_column{$effect} ]
                                )
                            )
                            && ( $transcript_effects
                                [ $vep_format_field_column{$effect} ] )
                          )
                        {

                            $vcf_info_csq_key{$effect}++;    #Increment
                        }
                    }
                }
            }
        }
    }
    close($FILEHANDLE);

    ##Check keys found in INFO field
    $log->info( 'Testing vcf INFO fields and presence in header: ' . $infile,
        "\n\n" );
    foreach my $key ( keys %vcf_info_key ) {

        ok(
            exists( $vcf_header{INFO}{$key} ),
            'Found both header and line field key for: '
              . $key
              . ' with key count: '
              . $vcf_info_key{$key}
        );
    }
    foreach my $key ( keys %vep_format_field_column ) {

        my @todo_keys = (
            'PolyPhen',  'APPRIS',
            'TSL',       'Existing_variation',
            'LoF_flags', 'MOTIF_NAME',
            'MOTIF_POS', 'HIGH_INF_POS',
            'MOTIF_SCORE_CHANGE', 'LoF_filter',
        );
        if ( any { $_ eq $key } @todo_keys ) {

            ## This will fail for Appris tcl etc which is only available in Grch38
          TODO: {
                local $TODO = 'Check VEP CSQ currently not produced';

                ok( $vcf_info_csq_key{$key},
                    'Found entry for CSQ field key for: ' . $key );
            }
        }
        else {

            ok( $vcf_info_csq_key{$key},
                    'Found entry for CSQ field key for: '
                  . $key
                  . ' with key count: '
                  . $vcf_info_csq_key{$key} );
        }
    }
    return;
}

sub read_range_file {

##read_range_file

##Function : Reads a file containg features to be annotated using range queries e.g. EnsemblGeneID.
##Returns  : ""
##Arguments: $vcfparser_data_href, $range_coulumns_ref, $range_file_key, $infile_path, $paddingRef, $select_feature_matching_column
##         : $vcfparser_data_href            => Range file hash {REF}
##         : $range_coulumns_ref             => Range columns to include {REF}
##         : $range_file_key                 => Range file key used to seperate range file(s) i.e., select and range
##         : $infile_path                    => Infile path
##         : $select_feature_matching_column => Column in the select file to match with vcf key annotation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vcfparser_data_href;
    my $range_coulumns_ref;
    my $range_file_key;
    my $infile_path;
    my $select_feature_matching_column;

    my $tmpl = {
        vcfparser_data_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$vcfparser_data_href
        },
        range_coulumns_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$range_coulumns_ref
        },
        range_file_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$range_file_key
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        select_feature_matching_column =>
          { strict_type => 1, store => \$select_feature_matching_column },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger('Test');

    my @headers;    #Save headers from rangeFile

    my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle
    open( $FILEHANDLE, '<', $infile_path )
      or $log->logdie( 'Cannot open ' . $infile_path . ':' . $!, "\n" );

    while (<$FILEHANDLE>) {

        chomp $_;                          #Remove newline

        if (m/^\s+$/) {                    # Avoid blank lines

            next;
        }
        if ( $_ =~ /^##/ ) {               #MetaData - Avoid

            next;
        }
        if ( $_ =~ /^#/ ) {                #Header/Comment

            @headers = split( /\t/, $_ );

            for (
                my $extract_columns_counter = 0 ;
                $extract_columns_counter < scalar(@$range_coulumns_ref) ;
                $extract_columns_counter++
              )
            {                              #Defines what scalar to store

                my $header_key_ref =
                  \$headers[ $$range_coulumns_ref[$extract_columns_counter] ];
                $vcfparser_data_href->{$$header_key_ref} =
                  $extract_columns_counter
                  ;    #Column position in supplied range input file
            }
            next;
        }
    }
    close($FILEHANDLE);
    $log->info(
        'Finished reading ' . $range_file_key . ' file: ' . $infile_path,
        "\n" );
    return;
}

sub parse_meta_data {

##parse_meta_data

##Function : Writes metadata to filehandle specified by order in meta_data_orders.
##Returns  : ""
##Arguments: $meta_data_href, $meta_data_string
##         : $meta_data_href   => Hash for meta_data {REF}
##         : $meta_data_string => The meta_data string from vcf header

    my ( $meta_data_href, $meta_data_string ) = @_;

    if ( $meta_data_string =~ /^##fileformat/ )
    {    #Catch fileformat as it has to be at the top of header

        push(
            @{ $meta_data_href->{fileformat}{fileformat} },
            $meta_data_string
        );    #Save metadata string
    }
    elsif ( $meta_data_string =~ /^##contig/ )
    {         #catch contigs to not sort them later

        push( @{ $meta_data_href->{contig}{contig} }, $meta_data_string )
          ;    #Save metadata string
    }
    elsif ( $meta_data_string =~ /^##(\w+)=(\S+)/ )
    {          #FILTER, FORMAT, INFO etc and more custom records

        push( @{ $meta_data_href->{$1}{$2} }, $meta_data_string )
          ;    #Save metadata string
    }
    else {     #All oddities

        push( @{ $meta_data_href->{other}{other} }, $meta_data_string )
          ;    #Save metadata string
    }
    return;
}

