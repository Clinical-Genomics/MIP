#!/usr/bin/env perl

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

##Collects MPS QC from MIP. Loads information on files to examine and values to extract from in YAML format and outputs exracted metrics in YAML format.
#Copyright 2013 Henrik Stranneheim

use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

## Third party module(s)
use YAML;

our $USAGE;

BEGIN {
    $USAGE =
        qq{qccollect.pl -si [sample_info.yaml] -r [regexp.yaml] -o [outfile]
               -si/--sample_info_file Sample info file (YAML format, Supply whole path, mandatory)
               -r/--regexp_file Regular expression file (YAML format, Supply whole path, mandatory)
               -o/--outfile The data file output (Supply whole path, defaults to "qcmetrics.yaml")
               -preg/--print_regexp Print the regexp used at CMMS switch (Defaults to "0" (=no))
               -prego/--print_regexp_outfile Regexp YAML outfile (Defaults to "qc_regexp.yaml")
               -h/--help Display this help message
               -v/--version Display version};
}
my ($sample_info_file, $regexp_file, $print_regexp);
my ($outfile, $print_regexp_outfile) = ("qcmetrics.yaml", "qc_regexp.yaml");
my (%qc_data, %evaluate_metric);
my %qc_header; #Save header(s) in each outfile
my %qc_program_data; #Save data in each outfile

my $qccollect_version = "2.0.0";

GetOptions('si|sample_info_file:s' => \$sample_info_file,
	   'r|regexp_file:s' => \$regexp_file,
	   'o|outfile:s'  => \$outfile,
	   'preg|print_regexp:n' => \$print_regexp,
	   'prego|print_regexp_outfile:s' => \$print_regexp_outfile,
	   'h|help' => sub { say STDOUT $USAGE; exit;},  #Display help text
	   'v|version' => sub { say STDOUT "\nqccollect.pl ".$qccollect_version, "\n"; exit;},  #Display version number
    )  or help({USAGE => $USAGE,
		exit_code => 1,
	       });

if ($print_regexp) {

    ## Write default regexp to YAML
    regexp_to_yaml({print_regexp_outfile => $print_regexp_outfile,
		   });
    print STDOUT "Wrote regexp YAML file to: ".$print_regexp_outfile, "\n";
    exit;
}

if (! $sample_info_file) {

    print STDERR "\n";
    print STDOUT $USAGE, "\n";
    print STDERR "Must supply a '-sample_info_file' (supply whole path)", "\n\n";
    exit;
}
if (! $regexp_file) {

    print STDERR "\n";
    print STDOUT $USAGE, "\n";
    print STDERR "Must supply a '-regexp_file' (supply whole path)", "\n\n";
    exit;
}

####MAIN

## Loads a YAML file into an arbitrary hash and returns it
my %sample_info = load_yaml({yaml_file => $sample_info_file,
			    });

## Loads a YAML file into an arbitrary hash and returns it
my %regexp = load_yaml({yaml_file => $regexp_file,
		       });


## Extracts all qcdata on sample_id level using information in %sample_info and %regexp
sample_qc({sample_info_href => \%sample_info,
	   regexp_href => \%regexp,
	   qc_data_href => \%qc_data,
	   qc_header_href => \%qc_header,
	   qc_program_data_href => \%qc_program_data,
	  });


## Extracts all qcdata on family level using information in %sample_info_file and %regexp
family_qc({sample_info_href => \%sample_info,
	   regexp_href => \%regexp,
	   qc_data_href => \%qc_data,
	   qc_header_href => \%qc_header,
	   qc_program_data_href => \%qc_program_data,
	  });

##Add qcCollect version to qc_data yaml file
$qc_data{program}{qccollect}{version} = $qccollect_version;
$qc_data{program}{qccollect}{regexp_file} = $regexp_file;

foreach my $sample_id (keys %{ $sample_info{sample} }) {

    ## Defines programs, etrics and thresholds to evaluate
    define_evaluate_metric({sample_info_href => \%sample_info,
			    sample_id => $sample_id,
			   });
}


## Evaluate the metrics
evaluate_qc_parameters({qc_data_href => \%qc_data,
			evaluate_metric_href => \%evaluate_metric,
		       });

## Writes a YAML hash to file
write_yaml({yaml_href => \%qc_data,
	    yaml_file_path_ref => \$outfile,
	   });

####SubRoutines

sub family_qc {

##family_qc

##Function : Extracts all qcdata on family level using information in %sample_info_file and %regexp
##Returns  : ""
##Arguments: $sample_info_href, $regexp_href, $qc_data_href, $qc_header_href, $qc_program_data_href
##         : $sample_info_href     => Info on samples and family hash {REF}
##         : $regexp_href          => RegExp hash {REF}
##         : $qc_data_href         => QCData hash {REF}
##         : $qc_header_href       => Save header(s) in each outfile {REF}
##         : $qc_program_data_href => Hash to save data in each outfile {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $regexp_href;
    my $qc_data_href;
    my $qc_header_href;
    my $qc_program_data_href;

    my $tmpl = {
	sample_info_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sample_info_href},
	regexp_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$regexp_href},
	qc_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_data_href},
	qc_header_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_header_href},
	qc_program_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_program_data_href},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    for my $program ( keys %{ $sample_info_href->{program} } ) { #For every program

	my $outdirectory;
	my $outfile;

	if ($sample_info_href->{program}{$program}{version} ) {

	    $qc_data_href->{program}{$program}{version} = $sample_info_href->{program}{$program}{version}; #Add version to qc_data
	}
	if ($sample_info_href->{program}{$program}{outdirectory} ) {

	    $outdirectory = $sample_info_href->{program}{$program}{outdirectory}; #Extract OutDirectory
	}
	if ($sample_info_href->{program}{$program}{outfile} ) {

	    $outfile = $sample_info_href->{program}{$program}{outfile}; #Extract OutFile
	}

	## Parses the RegExpHash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
	parse_regexp_hash_and_collect({regexp_href => $regexp_href,
				       qc_program_data_href => $qc_program_data_href,
				       qc_header_href => $qc_header_href,
				       program => $program,
				       outdirectory => $outdirectory,
				       outfile => $outfile,
				      });

	## Add extracted information to qc_data
	add_to_qc_data({sample_info_href => $sample_info_href,
			regexp_href => $regexp_href,
			qc_data_href => $qc_data_href,
			qc_header_href => $qc_header_href,
			qc_program_data_href => $qc_program_data_href,
			program => $program,
		       });
    }
}

sub sample_qc {

##sample_qc

##Function : Collects all sample qc in files defined by sample_info_file and regular expressions defined by regexp.
##Returns  : ""
##Arguments: $sample_info_href, $regexp_href, $qc_data_href, $qc_header_href, $qc_program_data_href
##         : $sample_info_href     => Info on samples and family hash {REF}
##         : $regexp_href          => RegExp hash {REF}
##         : $qc_data_href         => QCData hash {REF}
##         : $qc_header_href       => Save header(s) in each outfile {REF}
##         : $qc_program_data_href => Hash to save data in each outfile {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $regexp_href;
    my $qc_data_href;
    my $qc_header_href;
    my $qc_program_data_href;

    my $tmpl = {
	sample_info_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sample_info_href},
	regexp_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$regexp_href},
	qc_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_data_href},
	qc_header_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_header_href},
	qc_program_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_program_data_href},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

  SAMPLE_ID:
    for my $sample_id ( keys %{ $sample_info_href->{sample} } ) { #For every sample id

      PROGRAM:
	for my $program ( keys %{ $sample_info_href->{sample}{$sample_id}{program} } ) { #For every program

	  INFILE:
	    for my $infile ( keys %{ $sample_info_href->{sample}{$sample_id}{program}{$program} } ) { #For every infile

		my $outdirectory = $sample_info_href->{sample}{$sample_id}{program}{$program}{$infile}{outdirectory};
		my $outfile = $sample_info_href->{sample}{$sample_id}{program}{$program}{$infile}{outfile};

		## Parses the RegExpHash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
		parse_regexp_hash_and_collect({regexp_href => $regexp_href,
					       qc_program_data_href => $qc_program_data_href,
					       qc_header_href => $qc_header_href,
					       program => $program,
					       outdirectory => $outdirectory,
					       outfile => $outfile,
					      });

		## Add extracted information to qc_data
		add_to_qc_data({sample_info_href => $sample_info_href,
				regexp_href => $regexp_href,
				qc_data_href => $qc_data_href,
				qc_header_href => $qc_header_href,
				qc_program_data_href => $qc_program_data_href,
				sample_id => $sample_id,
				program => $program,
				infile => $infile,
			       });
	    }
	}
    }
}

sub parse_regexp_hash_and_collect {

##parse_regexp_hash_and_collect

##Function  : Parses the regexp hash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
##Returns   : ""
##Arguments : $regexp_href, qc_program_data_href, $qc_header_href, $program, $outdirectory, $outfile
##          : $regexp_href          => Regexp hash {REF}
##          : $qc_header_href       => Save header(s) in each outfile {REF}
##          : $qc_program_data_href => Hash to save data in each outfile {REF}
##          : $program              => The program to examine
##          : $outdirectory         => Programs out directory
##          : $outfile              => Programs out file containing parameter to evaluate

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regexp_href;
    my $qc_header_href;
    my $qc_program_data_href;
    my $program;
    my $outdirectory;
    my $outfile;

    my $tmpl = {
	regexp_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$regexp_href},
	qc_header_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_header_href},
	qc_program_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_program_data_href},
	program => { required => 1, defined => 1, strict_type => 1, store => \$program},
	outdirectory => { strict_type => 1, store => \$outdirectory},
	outfile => { strict_type => 1, store => \$outfile},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $regexp; #Holds the current regexp
    my @separators = ('\s+','!',','); #Covers both whitespace and tab. Add other separators if required

    for my $regexp_key ( keys %{ $regexp_href->{$program} } ) { #Find the actual regular expression(s) for each program that is used

	if ($regexp_key =~/^header|header$/i) { #Detect if the outfile contains paragrafs/header info in the outfile i.e. data is formated as a paragraf with header(s) and line(s). "regexp_key" should either start with or end with "header". This section extracts the header line(s) for the entire outdata file. Necessary to assign correct data entry to header entry later (headers and data are saved in seperate hashes).

	    ##Format outfile: Paragraf section
	  PARAGRAPH:
	    for my $regexp_header_key ( keys %{ $regexp_href->{$program}{$regexp_key} } ) { #Paragraf header

		$regexp = $regexp_href->{$program}{$regexp_key}{$regexp_header_key}; #The regular expression used to collect paragraf header info

	      SEPARATORS:
		for (my $separator_element_counter=0;$separator_element_counter<scalar(@separators);$separator_element_counter++) { #Loop through possible separators to seperate any eventual header elements

		    if ($regexp_header_key =~/^header|header$/i) { #Detect if the regexp key is a paragraf header and not paragraf data (header line and data line(s))

			@{ ${$qc_header_href}{$program}{$regexp_key}{$regexp_header_key} } = split(/$separators[$separator_element_counter]/, `$regexp $outdirectory/$outfile`); #Collect paragraf header

			if ( defined(${$qc_header_href}{$program}{$regexp_key}{$regexp_header_key})) { #Then split should have been successful

			    last; #Found correct separator - do not continue
			}
		    }
		    else { #For paragraf data line(s)

			@{ $qc_program_data_href->{$program}{$regexp_key}{$regexp_header_key} } = split(/$separators[$separator_element_counter]/, `$regexp $outdirectory/$outfile`); #Collect paragraf data

			if ( defined($qc_program_data_href->{$program}{$regexp_key}{$regexp_header_key}[1])) { #Then split should have been successful
			    last; #Found correct separator - do not continue
			}
		    }
		}
	    }
	}
	else { #For info contained in Entry --> Value i.e. same line.

	    $regexp = $regexp_href->{$program}{$regexp_key}; #The regular expression used to collect info

	    for (my $separator_element_counter=0;$separator_element_counter<scalar(@separators);$separator_element_counter++) { #Loop through possible separators

		@{ $qc_program_data_href->{$program}{$regexp_key} } = split(/$separators[$separator_element_counter]/, `$regexp $outdirectory/$outfile`); #Collect data. Use regexp_key as element header

		if ( defined($qc_program_data_href->{$program}{$regexp_key}[1])) { #Then split should have been successful

		    last; #Found correct separator do not continue
		}
	    }
	}
    }
}


sub add_to_qc_data {

##add_to_qc_data

##Function  : Add to qc_data hash to enable write to yaml format
##Returns   : ""
##Arguments : $sample_info_href, $regexp_href, $qc_data_href, $qc_header_href, $qc_program_data_href, $sample_id, $program, $inFile
##          : $sample_info_href     => Info on samples and family hash {REF}
##          : $regexp_href          => RegExp hash {REF}
##          : $qc_data_href         => QCData hash {REF}
##          : $qc_header_href       => Save header(s) in each outfile {REF}
##          : $qc_program_data_href => Hash to save data in each outfile {REF}
##          : $sample_id            => SampleID
##          : $program              => The program to examine
##          : $inFile               => infile to program

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $regexp_href;
    my $qc_data_href;
    my $qc_header_href;
    my $qc_program_data_href;
    my $sample_id;
    my $program;
    my $infile;

    my $tmpl = {
	sample_info_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sample_info_href},
	regexp_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$regexp_href},
	qc_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_data_href},
	qc_header_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_header_href},
	qc_program_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_program_data_href},
	sample_id => { strict_type => 1, store => \$sample_id},
	program => { required => 1, defined => 1, strict_type => 1, store => \$program},
	infile => { strict_type => 1, store => \$infile},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

  REGEXP:
    for my $regexp_key ( keys %{ $regexp_href->{$program} } ) { #All regexp per program

	if ($regexp_key !~/^header|header$/i) { #For info contained in entry --> Value i.e. same line

	    if (scalar(@{ $qc_program_data_href->{$program}{$regexp_key} }) == 1) { #Enable seperation of writing array or key-->value in qc_data

		if ( ($sample_id) && ($infile) ) {

		    $qc_data_href->{sample}{$sample_id}{$infile}{$program}{$regexp_key} = $qc_program_data_href->{$program}{$regexp_key}[0]; #key-->value for sample_id
		}
		else {  #Family level

		    $qc_data_href->{program}{$program}{$regexp_key} = $qc_program_data_href->{$program}{$regexp_key}[0]; #key-->value for familyID
		}
		if ($program eq "chanjo_sexcheck") {#Check gender for sample_id

		    my $chanjo_sexcheck = @{ $qc_program_data_href->{$program}{$regexp_key} }[0]; #Array_ref

		    ## Check that assumed gender is supported by coverage on chrX and chrY
		    gender_check({sample_info_href => $sample_info_href,
				  qc_data_href => $qc_data_href,
				  sample_id_ref => \$sample_id,
				  infile_ref => \$infile,
				  chanjo_sexcheck_gender_ref => \$chanjo_sexcheck,
				 });
		}
	    }
	    else { #Write array to qc_data

		for (my $regexp_key_counter=0;$regexp_key_counter<scalar(@{ $qc_program_data_href->{$program}{$regexp_key} });$regexp_key_counter++ ) {

		    if ( ($sample_id) && ($infile) ) {

			$qc_data_href->{sample}{$sample_id}{$infile}{$program}{$regexp_key}[$regexp_key_counter] = $qc_program_data_href->{$program}{$regexp_key}[$regexp_key_counter];

		    }
		    else {

			$qc_data_href->{program}{$program}{$regexp_key}[$regexp_key_counter] = $qc_program_data_href->{$program}{$regexp_key}[$regexp_key_counter];
		    }
		    if ($program eq "plink_sexcheck") {#Check gender for sample_id

			my @sexchecks = split(":", @{ $qc_program_data_href->{$program}{$regexp_key} }[$regexp_key_counter]); #ArrayRef

			## Check that assumed gender is supported by variants on chrX and chrY
			plink_gender_check({sample_info_href => $sample_info_href,
					    qc_data_href => $qc_data_href,
					    sample_id_ref => \$sexchecks[0],
					    plink_sexcheck_gender_ref => \$sexchecks[1],
					   });
		    }
		}
		if (defined($qc_data_href->{program}{relation_check}{sample_relation_check}) && (defined($qc_data_href->{program}{pedigree_check}{sample_order}) ) ) {

		    relation_check({sample_info_href => $sample_info_href,
				    qc_data_href => $qc_data_href,
				    relationship_values_ref => \@{ $qc_data_href->{program}{relation_check}{sample_relation_check} },
				    sample_orders_ref => \@{ $qc_data_href->{program}{pedigree_check}{sample_order} },
				   });
		}
	    }
	}
	else { #Paragraf data i.e. header and subsequent data lines

	  HEADER_INFO:
	    for my $regexp_header_key ( keys %{ ${$qc_header_href}{$program}{$regexp_key} }) { #Find header info

	      PARAGRAPH_KEYS:
		for my $regexp_key_header ( keys %{ $regexp_href->{$program}{$regexp_key} } ) { #All paragraf keys (header and data line(s))

		    if ($regexp_key_header !~/^header|header$/i) { #Detect if the regexp id for headers and not data.

			for (my $qc_headers_counter=0;$qc_headers_counter<scalar( @{ ${$qc_header_href}{$program}{$regexp_key}{$regexp_header_key} } );$qc_headers_counter++) { #For all collected headers

			    if ( ($sample_id) && ($infile)) {

				$qc_data_href->{sample}{$sample_id}{$infile}{$program}{$regexp_header_key}{$regexp_key_header}{ ${$qc_header_href}{$program}{$regexp_key}{$regexp_header_key}[$qc_headers_counter] } = $qc_program_data_href->{$program}{$regexp_key}{$regexp_key_header}[$qc_headers_counter]; #Add to qc_data using header element[X] --> data[X] to correctly position elements in qc_data hash
                            }
                            else {

				$qc_data_href->{$program}{$regexp_header_key}{$regexp_key_header}{ ${$qc_header_href}{$program}{$regexp_key}{$regexp_header_key}[$qc_headers_counter] } = $qc_program_data_href->{$program}{$regexp_key}{$regexp_key_header}[$qc_headers_counter]; #Add to qc_data using header element[X] --> data[X] to correctly position elements in qc_data hash

                            }
                        }
		    }
		}
	    }
	}
    }
}

sub define_evaluate_metric {

##define_evaluate_metric

##Function  : Sets programs and program metrics and thresholds to be evaluated
##Returns   : ""
##Arguments : $sample_info_href, $sample_id
##          : $sample_info_href => Info on samples and family hash {REF}
##          : $sample_id        => SampleID

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $sample_id;

    my $tmpl = {
	sample_info_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sample_info_href},
	sample_id => { required => 1, defined => 1, strict_type => 1, store => \$sample_id},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    $evaluate_metric{$sample_id}{mosaik_aligner}{total_aligned}{threshold} = 95;
    $evaluate_metric{$sample_id}{mosaik_aligner}{uniquely_aligned_mates}{threshold} = 90;
    $evaluate_metric{$sample_id}{bamstats}{percentage_mapped_reads}{threshold} = 0.95;
    $evaluate_metric{$sample_id}{calculatehsmetrics}{PCT_TARGET_BASES_10X}{threshold} = 0.95;
    $evaluate_metric{$sample_id}{collectmultiplemetrics}{PCT_PF_READS_ALIGNED}{threshold} = 0.95;
    $evaluate_metric{$sample_id}{calculatehsmetrics}{PCT_ADAPTER}{threshold} = 0.0001;

    if (exists($sample_info_href->{sample}{$sample_id}{expected_coverage})) {

	$evaluate_metric{$sample_id}{calculatehsmetrics}{MEAN_TARGET_COVERAGE}{threshold} = $sample_info_href->{sample}{$sample_id}{expected_coverage};
    }

    $evaluate_metric{variant_integrity_mendel}{fraction_of_errors}{gt} = 0.06;
    $evaluate_metric{variant_integrity_father}{fraction_of_common_variants}{lt} = 0.55;

    if ($sample_info{sample}{$sample_id}{analysis_type} eq "wes") {

	$evaluate_metric{$sample_id}{calculatehsmetrics}{PCT_TARGET_BASES_30X}{threshold} = 0.90;
    }
}
sub evaluate_qc_parameters {

##evaluate_qc_parameters

##Function  : Evaluate parameters to detect parameters falling below threshold
##Returns   : ""
##Arguments : $qc_data_href, $evaluate_metric_href
##          : $qc_data_href         => QCData hash {REF}
##          : $evaluate_metric_href => HAsh for metrics to evaluate

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $evaluate_metric_href;

    my $tmpl = {
	qc_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_data_href},
	evaluate_metric_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$evaluate_metric_href},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $status;

  PROGRAM:
    for my $program ( keys %{$qc_data_href->{program}} ) {

	if (defined($evaluate_metric_href->{$program})) { #Program to be evaluated

	  METRIC:
	    for my $metric ( keys %{$qc_data_href->{program}{$program}} ) {

		if (defined($evaluate_metric_href->{$program}{$metric})) { #metric to be evaluated

		    if ($evaluate_metric_href->{$program}{$metric}{gt}) {

			if ($qc_data_href->{program}{$program}{$metric} > $evaluate_metric_href->{$program}{$metric}{gt}) { #Determine status - if greater than add to hash. otherwise PASS and do not include

			    $status = "FAILED:".$program."_".$metric.":".$qc_data_href->{program}{$program}{$metric};
			    push(@{$qc_data_href->{evaluation}{$program}}, $status);
			}
		    }
		    if ($evaluate_metric_href->{$program}{$metric}{lt}) {

			if ($qc_data_href->{program}{$program}{$metric} < $evaluate_metric_href->{$program}{$metric}{lt}) { #Determine status - if lower than add to hash. otherwise PASS and do not include

			    $status = "FAILED:".$program."_".$metric.":".$qc_data_href->{program}{$program}{$metric};
			    push(@{$qc_data_href->{evaluation}{$program}}, $status);
			}
		    }
		    last;
		}
	    }
	}
    }

    ## Sample level evaluation
  SAMPLE_ID:
    for my $sample_id ( keys %{$qc_data_href->{sample}} ) {

      INFILE:
	for my $infile ( keys %{$qc_data_href->{sample}{$sample_id}} ) {

	    if ($infile =~/relation_check/) { #Special case

		if ($qc_data_href->{sample}{$sample_id}{$infile} ne "PASS") {

		    $status = "Status:".$infile.":".$qc_data_href->{sample}{$sample_id}{$infile};
		    push(@{$qc_data_href->{evaluation}{$infile}}, $status); #Add to QC data at family level
		}
		next;
	    }
	    if ($infile =~/evaluation/) { #Special case

		next;
	    }
	  PROGRAM:
	    for my $program ( keys %{$qc_data_href->{sample}{$sample_id}{$infile}} ) {

		if (defined($evaluate_metric_href->{$sample_id}{$program})) { #Program to be evaluated

		  METRIC:
		    for my $metric ( keys %{$evaluate_metric_href->{$sample_id}{$program}}) { #Metric to be evaluated

			if (defined($qc_data_href->{sample}{$sample_id}{$infile}{$program}{$metric})) {

			    if ($qc_data_href->{sample}{$sample_id}{$infile}{$program}{$metric} < $evaluate_metric_href->{$sample_id}{$program}{$metric}{threshold}) { #Determine status - if below add to hash. otherwise PASS and do not include

				$status = "FAILED:".$sample_id."_".$program."_".$metric.":".$qc_data_href->{sample}{$sample_id}{$infile}{$program}{$metric};
				push(@{$qc_data_href->{evaluation}{$program}}, $status);
			    }
			    last;
			}
			else {

			    for my $key ( keys %{$qc_data_href->{sample}{$sample_id}{$infile}{$program}} ) {

				if ($key eq "header") {

				    for my $data_header ( keys %{$qc_data_href->{sample}{$sample_id}{$infile}{$program}{$key}} ) {

					if (defined($qc_data_href->{sample}{$sample_id}{$infile}{$program}{$key}{$data_header}{$metric})) {

					    if ($qc_data_href->{sample}{$sample_id}{$infile}{$program}{$key}{$data_header}{$metric} < $evaluate_metric_href->{$sample_id}{$program}{$metric}{threshold}) { #Determine status - if below add to hash. otherwise PASS and do not include

						$status = "FAILED:".$sample_id."_".$program."_".$metric.":".$qc_data_href->{sample}{$sample_id}{$infile}{$program}{$key}{$data_header}{$metric};
						push(@{$qc_data_href->{evaluation}{$program}}, $status);
					    }
					    next; #Metric go to next section
					}
				    }
				    last; #Metric found no need to continue
				}
			    }
			}
		    }
		}
	    }
	}
    }
}


sub relation_check {

##relation_check

##Function : Uses the .mibs file produced by PLINK to test if family members are indeed related.
##Returns  : ""
##Arguments: $sample_info_href, $qc_data_href, $relationship_values_ref, $sample_orders_ref
##         : $sample_info_href        => Info on samples and family hash {REF}
##         : $qc_data_href            => QCData hash {REF}
##         : $sample_orders_ref       => The sample order so that correct estimation can be connected to the correct sample_ids {REF}
##         : $relationship_values_ref => All relationship estimations {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $qc_data_href;
    my $relationship_values_ref;
    my $sample_orders_ref;

    my $tmpl = {
	sample_info_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sample_info_href},
	qc_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_data_href},
	relationship_values_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$relationship_values_ref},
	sample_orders_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$sample_orders_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my %family; #Stores family relations and pairwise comparisons family{$sample_id}{$sample_id}["column"] -> [pairwise]
    my $sample_id_counter = 0;
    my $incorrect_relation = 0;
    my @pairwise_comparisons;
    my @relationship_values = @$relationship_values_ref;  #Copy array to avoid removing actual values in later splice

    ## Splice all relationship extimations from regexp into pairwise comparisons calculated for each sample_id
    for (my $realtionship_counter=0;$realtionship_counter<scalar(@relationship_values);$realtionship_counter++) {

	my @pairwise_comparisons = splice(@relationship_values, 0, scalar(@$sample_orders_ref)); #Splices array into each sample_ids line

	for (my $column=0;$column<scalar(@$sample_orders_ref);$column++) { #All columns in .mibs file

	    push ( @{ $family{ $sample_orders_ref->[$sample_id_counter] }{ $sample_orders_ref->[$column] } }, $pairwise_comparisons[$column]); #Store sample_id, family membersID (including self) and each pairwise comparison. Uses array for to accomodate sibling info.
	}
	$sample_id_counter++;
    }
    my $father_id = "YYY"; #father_id for the family
    my $mother_id = "XXX"; #mother_id for the family

    ## Collect father and mother id
    for my $sample_id ( keys %family ) { #For all sample_ids

	## Currently only 1 father or Mother per pedigree is supported

	if ($sample_info_href->{sample}{$sample_id}{father} ne 0) { #Save father_id if not 0

	    $father_id = $sample_info_href->{sample}{$sample_id}{father};
	}
	if ($sample_info_href->{sample}{$sample_id}{mother} ne 0) { #Save mother_id if not 0

	    $mother_id = $sample_info_href->{sample}{$sample_id}{mother};
	}
    }

  SAMPLE_ID:
    for my $sample_id ( keys %family ) { #For all sample_ids

      MEMBER:
	for my $members ( keys %{ $family{$sample_id} } ) { #For every relation within family (mother/father/child)

	  RELATIVES:
	    for (my $members_count=0;$members_count<scalar( @{ $family{$sample_id}{$members} } );$members_count++) { #@ Necessary for siblings

		if ($family{$sample_id}{$members}[$members_count] == 1 ) { #Should only hit self

		    if ( $sample_id eq  $members) {

			#print "Self: ".$sample_id,"\t", $members, "\t", $family{$sample_id}{$members}[$members_count], "\n";
		    }
		    else {

			$incorrect_relation++;
			$qc_data_href->{sample}{$sample_id}{relation_check} = "FAIL: Duplicated sample?;";
			#print  "Incorrect should be self: ".$sample_id,"\t", $members, "\t", $family{$sample_id}{$members}[$members_count], "\n";
		    }
		}
		elsif ($family{$sample_id}{$members}[$members_count] >= 0.70 ) { #Should include parent to child and child to siblings unless inbreed parents

		    if ( ( ($sample_id ne $father_id) && ($sample_id ne $mother_id) )
			 || ( ($members ne $father_id) && ($members ne $mother_id) ) ) { #Correct
			#print "Parent-to-child or child-to-child: ".$sample_id,"\t", $members, "\t", $family{$sample_id}{$members}[$members_count], "\n";
		    }
		    else {

			$incorrect_relation++;
			$qc_data_href->{sample}{$sample_id}{relation_check} = "FAIL: Parents related?;";
			#print "Incorrect: ".$sample_id,"\t", $members, "\t", $family{$sample_id}{$members}[$members_count], "\n";
		    }
		}
		elsif ($family{$sample_id}{$members}[$members_count] < 0.70 ) { #Parents unless inbreed

		    if ( ($sample_id eq $father_id) && ($members eq $mother_id) ) {

			#print "Parents: ".$sample_id,"\t", $members, "\t", $family{$sample_id}{$members}[$members_count], "\n";
		    }
		    elsif ( ($sample_id eq $mother_id) && ($members eq $father_id) ) {

			#print "Parents: ".$sample_id,"\t", $members, "\t", $family{$sample_id}{$members}[$members_count], "\n";
		    }
		    else {

			$incorrect_relation++;
			$qc_data_href->{sample}{$sample_id}{relation_check} = "FAIL:".$sample_id." not related to ".$members.";";
			#print "Incorrect: ".$sample_id,"\t", $members, "\t", $family{$sample_id}{$members}[$members_count], "\n";
		    }
		}
	    }
	}
	if ($incorrect_relation == 0) {

	    $qc_data_href->{sample}{$sample_id}{relation_check} = "PASS";
	}
    }
    return;
}

sub gender_check {

##Function : Checks that the gender predicted by chanjo_sexcheck is confirmed in the pedigee for the sample
##Returns  : ""
##Arguments: $sample_info_href, $qc_data_href, sample_id_ref, $infile_ref, $chanjo_sexcheck_gender_ref
##         : $sample_info_href           => Info on samples and family hash {REF}
##         : $qc_data_href               => QCData hash {REF}
##         : $sample_id_ref              => SampleID {REF}
##         : $infile_ref                 => Infile {REF}
##         : $chanjo_sexcheck_gender_ref => Chanjo calculated gender {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $qc_data_href;
    my $sample_id_ref;
    my $infile_ref;
    my $chanjo_sexcheck_gender_ref;

    my $tmpl = {
	sample_info_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sample_info_href},
	qc_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_data_href},
	sample_id_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$sample_id_ref},
	infile_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$infile_ref},
	chanjo_sexcheck_gender_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$chanjo_sexcheck_gender_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $sample_id_sex_ref = \$sample_info_href->{sample}{$$sample_id_ref}{sex}; #Alias

    if ( ($$chanjo_sexcheck_gender_ref eq "female") && ($$sample_id_sex_ref =~/2|female/) ) { #Female

	$qc_data_href->{sample}{$$sample_id_ref}{$$infile_ref}{gender_check} = "PASS";
    }
    elsif ( ($$chanjo_sexcheck_gender_ref eq "male") && ($$sample_id_sex_ref =~/1|^male/) ) { #Male

	$qc_data_href->{sample}{$$sample_id_ref}{$$infile_ref}{gender_check} = "PASS";
    }
    else {

	$qc_data_href->{sample}{$$sample_id_ref}{$$infile_ref}{gender_check} = "FAIL";
    }
    return;
}


sub plink_gender_check {

##plink_gender_check

##Function : Checks that the gender predicted by Plink sexcheck is confirmed in the pedigee for the sample
##Returns  : ""
##Arguments: $sample_info_href, $qc_data_href, $sample_id_ref, $plink_sexcheck_gender_ref
##         : $sample_info_href          => Info on samples and family hash {REF}
##         : $qc_data_href              => QCData hash {REF}
##         : $sample_id_ref             => SampleID {REF}
##         : $plink_sexcheck_gender_ref => Plink calculated gender {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $qc_data_href;
    my $sample_id_ref;
    my $plink_sexcheck_gender_ref;

    my $tmpl = {
	sample_info_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sample_info_href},
	qc_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qc_data_href},
	sample_id_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$sample_id_ref},
	plink_sexcheck_gender_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$plink_sexcheck_gender_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $sample_id_sex_ref = \$sample_info_href->{sample}{$$sample_id_ref}{sex}; #Alias

    if ( ($$plink_sexcheck_gender_ref eq "2") && ($$sample_id_sex_ref =~/2|female/) ) { #Female

	push(@{$qc_data_href->{program}{plink_gender_check}}, $$sample_id_ref.":PASS");
    }
    elsif ( ($$plink_sexcheck_gender_ref eq "1") && ($$sample_id_sex_ref =~/1|^male/) ) { #Male

	push(@{$qc_data_href->{program}{plink_gender_check}}, $$sample_id_ref.":PASS");
    }
    else {

	push(@{$qc_data_href->{program}{plink_gender_check}}, $$sample_id_ref.":FAIL");
    }
    return;
}

sub write_yaml {

##write_yaml

##Function : Writes a YAML hash to file
##Returns  : ""
##Arguments: $yaml_href, $yaml_file_path_ref
##         : $yaml_href          => The hash to dump {REF}
##         : $yaml_file_path_ref => The yaml file to write to {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $yaml_href;
    my $yaml_file_path_ref;

    my $tmpl = {
	yaml_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$yaml_href},
	yaml_file_path_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$yaml_file_path_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    open (my $YAML, ">", $$yaml_file_path_ref) or die "can't open ".$$yaml_file_path_ref.": $!\n";
    say $YAML Dump( $yaml_href );
    close($YAML);
}

sub load_yaml {

##load_yaml

##Function : Loads a YAML file into an arbitrary hash and returns it.
##Returns  : %yaml_hash
##Arguments: $yaml_file
##         : $yaml_file => The yaml file to load

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $yaml_file;

    my $tmpl = {
	yaml_file => { required => 1, defined => 1, strict_type => 1, store => \$yaml_file},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my %yaml;

    open (my $YAML, "<", $yaml_file) or die "can't open ".$yaml_file.":".$!, "\n";

    %yaml = %{ YAML::LoadFile($yaml_file) };  #Load hashreference as hash

    close($YAML);

    return %yaml;
}


sub regexp_to_yaml {

##regexp_to_yaml

##Function : Write default regexp to YAML
##Returns  : ""
##Arguments: $print_regexp_outfile
##         : $print_regexp_outfile => File to print regexp to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $print_regexp_outfile;

    my $tmpl = {
	print_regexp_outfile => { required => 1, defined => 1, strict_type => 1, store => \$print_regexp_outfile},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my %regexp;

    #Add to %regexp to enable print in YAML
    $regexp{fastqc}{version} = q?perl -nae' if ($_=~/##FastQC\\s+(\\S+)/) {print $1;last;}' ?; #Collect FastQC version

    $regexp{fastqc}{encoding} = q?perl -nae' if ($_=~/Encoding\s+(\S+\s\S+\s\S+\s\S+|\S+\s\S+)/) { my $encoding = $1;$encoding=~s/\s/\_/g; print $encoding;last;}' ?; #Collect Encoding

    $regexp{fastqc}{sequence_length} = q?perl -nae' if ($_=~/Sequence length\s(\d+)/) {print $1;last;}' ?; #Collect Sequence length

    $regexp{fastqc}{total_number_of_reads} = q?perl -nae' if ($_=~/Total Sequences\s(\d+)/) {print $1;last;}' ?; #Collect Total sequences

    $regexp{fastqc}{gc} = q?perl -nae' if ($_=~/%GC\s(\d+)/) {print $1;last;}' ?; #Collect GC content

    $regexp{fastqc}{sequence_duplication} = q?perl -nae' if ($_=~/#Total Duplicate Percentage\s+(\d+.\d)/) {print $1;last;}' ?; #Collect Sequence duplication level

    $regexp{fastqc}{basic_statistics} = q?perl -nae' if ($_=~/>>Basic Statistics\s+(\S+)/) {print $1;last;}' ?; #Collect Basic Statistics

    $regexp{fastqc}{per_base_sequence_quality} = q?perl -nae' if ($_=~/>>Per base sequence quality\s+(\S+)/) {print $1;last;}' ?; #Collect Per base sequence quality

    $regexp{fastqc}{per_sequence_quality_scores} = q?perl -nae' if ($_=~/>>Per sequence quality scores\s+(\S+)/) {print $1;last;}' ?; #Collect Per sequence quality scores

    $regexp{fastqc}{per_base_sequence_content} = q?perl -nae' if ($_=~/>>Per base sequence content\s+(\S+)/) {print $1;last;}' ?; #Collect Per base sequence content

    $regexp{fastqc}{per_base_gc_content} = q?perl -nae' if ($_=~/>>Per base GC content\s+(\S+)/) {print $1;last;}' ?; #Collect Per base GC content

    $regexp{fastqc}{per_sequence_gc_content} = q?perl -nae' if ($_=~/>>Per sequence GC content\s+(\S+)/) {print $1;last;}' ?; #Collect Per sequence GC content

    $regexp{fastqc}{per_base_n_content} = q?perl -nae' if ($_=~/>>Per base N content\s+(\S+)/) {print $1;last;}' ?; #Collect Per base N content

    $regexp{fastqc}{sequence_duplication_levels} = q?perl -nae' if ($_=~/>>Sequence Duplication Levels\s+(\S+)/) {print $1;last;}' ?; #Collect Sequence Duplication Levels

    $regexp{fastqc}{overrepresented_sequences} = q?perl -nae' if ($_=~/>>Overrepresented sequences\s+(\S+)/) {print $1;last;}' ?; #Collect Overrepresented sequences

    $regexp{fastqc}{kmer_content} = q?perl -nae' if ($_=~/>>Kmer Content\s+(\S+)/) {print $1;last;}' ?; #Collect Kmer Content

    $regexp{mosaik_aligner}{version} = q?perl -nae' if ($_=~/(\d+\.\d+\.\d+)\s/) {print $1;last;}' ?; #Collect Mosaik Version

    $regexp{mosaik_aligner}{unaligned_mates} = q?perl -nae' if ($_=~/# unaligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Nr of unaligned mates

    $regexp{mosaik_aligner}{filtered_out} = q?perl -nae' if ($_=~/# filtered out\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Nr of filtered out reads

    $regexp{mosaik_aligner}{uniquely_aligned_mates} = q?perl -nae' if ($_=~/# uniquely aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Uniquely aligned mates

    $regexp{mosaik_aligner}{multiply_aligned_mates} = q?perl -nae' if ($_=~/# multiply aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Multiply aligned mates

    $regexp{mosaik_aligner}{total_aligned} = q?perl -nae' if ($_=~/total aligned:\s+\S+\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) {print $2;last;} elsif ($_=~/total aligned:\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) { print $2;last;}' ?; #Collect total aligned sequences

    $regexp{bamstats}{percentage_mapped_reads} = q?perl -nae 'if($_=~/percentage mapped reads:\s+(\S+)/) {print $1;last}' ?; #Collect % mapped reads from BAM alignment

    $regexp{bamstats}{raw_total_sequences} = q?perl -nae 'if($_=~/raw total sequences:\s+(\S+)/) {print $1;last}' ?; #Collect raw total sequences from BAM alignment

    $regexp{bamstats}{reads_mapped} = q?perl -nae 'if($_=~/reads mapped:\s+(\S+)/) {print $1;last}' ?; #Collect reads mapped from BAM alignment

    $regexp{chanjo_sexcheck}{gender} = q?perl -nae 'if( ($F[0]!~/^#/) && ($F[2] =~/\S+/) ) {print $F[2];}' ?;  #Collect gender from chanjo_sexcheck

    $regexp{pedigree_check}{sample_order} = q?perl -nae 'if ($_=~/^#CHROM/) {chomp $_; my @line = split(/\t/,$_); for (my $sample=9;$sample<scalar(@line);$sample++) { print $line[$sample], "\t";}last;}' ?; #Collect sample order from vcf file used to create ".ped", ".map" and hence ".mibs".

    $regexp{inbreeding_factor}{sample_inbreeding_factor}  = q?perl -nae 'my @inbreedingFactor; if ($. > 1) {my @temp = split(/\s/,$_);push(@inbreedingFactor, $F[0].":".$F[5]); print $inbreedingFactor[0], "\t"; }' ?;

    $regexp{plink_sexcheck}{sample_sexcheck}  = q?perl -nae 'my @sexCheckFactor; if ($. > 1) {my @temp = split(/\s+/,$_);push(@sexCheckFactor,$temp[2].":".$temp[4]); print $sexCheckFactor[0], "\t"; }' ?;

    $regexp{relation_check}{sample_relation_check}  = q?perl -nae 'print $_;' ?; #Note will return whole file

    $regexp{markduplicates}{fraction_duplicates} = q?perl -nae 'if($_=~/Fraction Duplicates\: (\S+)/) {print $1;}' ?; #Collect fraction duplicates

    $regexp{calculatehsmetrics}{header_info}{header} = q?perl -nae' if ($_ =~/^BAIT_SET/ ) {print $_;last;}' ?; #Note return whole line (header)

    $regexp{calculatehsmetrics}{header_info}{data} = q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?; #Note return whole line and only look at line 8, where the data action is

    $regexp{collectmultiplemetrics}{header_info}{header} = q?perl -nae' if ($_ =~/^CATEGORY/ ) {print $_;last;}' ?; #Note return whole line (header)

    $regexp{collectmultiplemetrics}{header_info}{first_of_pair} = q?perl -nae' if ($_ =~/^FIRST_OF_PAIR/ ) {print $_;last;}' ?; #Note return whole line (FIRST_OF_PAIR)

    $regexp{collectmultiplemetrics}{header_info}{second_of_pair} = q?perl -nae' if ($_ =~/^SECOND_OF_PAIR/ ) {print $_;last;}' ?; #Note return whole line (SECOND_OF_PAIR)

    $regexp{collectmultiplemetrics}{header_info}{pair} = q?perl -nae' if ($_ =~/^PAIR/ ) {print $_;last;}'  ?; #Note return whole line (PAIR)

    $regexp{collectmultiplemetricsinsertsize}{header_info}{header} = q?perl -nae' if ($_ =~/^MEDIAN_INSERT_SIZE/ ) {print $_;last;}' ?; #Note return whole line (header)

    $regexp{collectmultiplemetricsinsertsize}{header_info}{data} = q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?; #Note return whole line and only look at line 8, where the data action is

    $regexp{variantevalall}{comp_overlap_header}{comp_overlap_header} = q?perl -nae' if ($_ =~/^CompOverlap\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (header)

    $regexp{variantevalall}{comp_overlap_header}{comp_overlap_data_all} = q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/all/) && ($_ =~/none/)) {print $_;last;}' ?; #Note return whole line

    $regexp{variantevalall}{comp_overlap_header}{comp_overlap_data_known} = q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line

    $regexp{variantevalall}{comp_overlap_header}{comp_overlap_data_novel} = q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regexp{variantevalall}{count_variants_header}{count_variants_header} = q?perl -nae' if ($_ =~/^CountVariants\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (header)
    $regexp{variantevalall}{count_variants_header}{count_variants_data_all} = q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{count_variants_header}{count_variants_data_known} = q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{count_variants_header}{count_variants_data_novel} = q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regexp{variantevalall}{indel_summary_header}{indel_summary_header} = q?perl -nae' if ($_ =~/^IndelSummary\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (header)
    $regexp{variantevalall}{indel_summary_header}{indel_summary_data_all} = q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{indel_summary_header}{indel_summary_data_known} = q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{indel_summary_header}{indel_summary_data_novel} = q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regexp{variantevalall}{multiallelic_summary_header}{multiallelic_summary_header} = q?perl -nae' if ($_ =~/^MultiallelicSummary\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (header)
    $regexp{variantevalall}{multiallelic_summary_header}{multiallelic_summary_data_all} = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{multiallelic_summary_header}{multiallelic_summary_data_known} = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{multiallelic_summary_header}{multiallelic_summary_data_novel} = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regexp{variantevalall}{titv_variant_evaluator_header}{titv_variant_evaluator_header} = q?perl -nae' if ($_ =~/^TiTvVariantEvaluator\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (header)
    $regexp{variantevalall}{titv_variant_evaluator_header}{titv_variant_evaluator_data_all} = q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{titv_variant_evaluator_header}{titv_variant_evaluator_data_known} = q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{titv_variant_evaluator_header}{titv_variant_evaluator_data_novel} = q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regexp{variantevalall}{validation_report_header}{validation_report_header} = q?perl -nae' if ($_ =~/^ValidationReport\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (header)
    $regexp{variantevalall}{validation_report_header}{validation_report_data_all} = q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/all\s/) && ($_ =~/none\s/)) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{validation_report_header}{validation_report_data_known} = q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{validation_report_header}{validation_report_data_novel} = q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regexp{variantevalall}{variant_summary_header}{variant_summary_header} = q?perl -nae' if ($_ =~/^VariantSummary\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (header)
    $regexp{variantevalall}{variant_summary_header}{variant_summary_data_all} = q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{variant_summary_header}{variant_summary_data_known} = q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line
    $regexp{variantevalall}{variant_summary_header}{variant_summary_data_novel} = q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regexp{variantevalexome} = $regexp{variantevalall};

    $regexp{genmod}{version} = q?perl -nae 'if($_=~/##Software=<ID=genmod,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect Genmod version

    $regexp{snpeff}{version} = q?perl -nae 'if($_=~/##SnpSiftVersion=\"(.+),/) {my $ret=$1; $ret=~s/\s/_/g;print $ret;last;}' ?; #Collect SnpEff version

    $regexp{varianteffectpredictor}{version} = q?perl -nae 'if($_=~/##VEP=(\w+)/) {print $1;last;}' ?; #Collect varianteffectpredictor version

    $regexp{varianteffectpredictor}{cache} = q?perl -nae 'if($_=~/##VEP=\w+\s+cache=(\S+)/) {print $1;last;}' ?; #Collect varianteffectpredictor cache directory

    $regexp{varianteffectpredictor}{polyphen} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/polyphen=(\S+)/) {print $1;last;}' ?; #Collect varianteffectpredictor polyPhen version

    $regexp{varianteffectpredictor}{sift} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/sift=sift(\S+)/) {print $1;last;}' ?; #Collect varianteffectpredictor sift version

    $regexp{varianteffectpredictor}{gene_build} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/genebuild=(\S+)/) {print $1;last;}' ?; #Collect varianteffectpredictor geneBuild

    $regexp{varianteffectpredictor}{assembly} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/assembly=(\S+)/) {print $1;last;}' ?; #Collect varianteffectpredictor assembly

    $regexp{varianteffectpredictor}{hgmd_public} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/HGMD-PUBLIC=(\S+)/) {print $1;last;}' ?; #Collect varianteffectpredictor HGMD-PUBLIC version

    $regexp{varianteffectpredictor}{reg_build} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/regbuild=(\S+)/) {print $1;last;}' ?; #Collect varianteffectpredictor regbuild version

    $regexp{varianteffectpredictor}{gencode} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/gencode=\S+\s+(\d+)/) {print $1;last;}' ?; #Collect varianteffectpredictor gencode version

    $regexp{vcfparser}{version} = q?perl -nae 'if($_=~/##Software=<ID=vcfParser.pl,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect vcfparser version

    $regexp{bwa}{version} = q?perl -nae 'if($_=~/\[main\]\sVersion:\s(\S+)/) {print $1;last;}' ?; #Collect Bwa version

    $regexp{chanjo}{version} = q?perl -nae 'if($_=~/version\s(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect Chanjo version

    $regexp{vt}{version} = q?perl -nae 'if($_=~/decompose\sv(\S+)/) {print $1;last;}' ?; #Collect vt version

    $regexp{samtools}{version} = q?perl -nae 'if($_=~/samtoolsVersion=(\S+)/) {print $1;last;}' ?; #Collect Samtools version

    $regexp{bcftools}{version} = q?perl -nae 'if($_=~/bcftools_\w+Version=(\S+)/) {print $1;last;}' ?; #Collect Bcftools version

    $regexp{freebayes}{version} = q?perl -nae 'if($_=~/source=freeBayes\s(\S+)/) {print $1;last;}' ?; #Collect Freebayes version

    $regexp{delly}{version} = q?perl -nae 'if($_=~/SVMETHOD=EMBL\.DELLY(v\d+\.\d+\.\d+)/) {print $1;last }' ?; #Collect Delly version

    $regexp{manta}{version} = q?perl -nae 'if($_=~/GenerateSVCandidates\s+(\S+)/) {print $1;last}' ?; #Collect Manta version

    $regexp{sv_combinevariantcallsets}{vcfanno} = q?perl -nae 'if($_=~/vcfanno\sversion\s(\S+)/) {print $1;last;}' ?; #Collect SVVCFAnno version

    $regexp{sv_varianteffectpredictor}{version} = q?perl -nae 'if($_=~/##VEP=(\w+)/) {print $1;last;}' ?; #Collect sv_varianteffectpredictor version

    $regexp{sv_varianteffectpredictor}{cache} = q?perl -nae 'if($_=~/##VEP=\w+\s+cache=(\S+)/) {print $1;last;}' ?; #Collect sv_varianteffectpredictor cache directory

    $regexp{sv_varianteffectpredictor}{polyphen} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/polyphen=(\S+)/) {print $1;last;}' ?; #Collect sv_varianteffectpredictor polyPhen version

    $regexp{sv_varianteffectpredictor}{sift} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/sift=sift(\S+)/) {print $1;last;}' ?; #Collect sv_varianteffectpredictor sift version

    $regexp{sv_varianteffectpredictor}{gene_build} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/genebuild=(\S+)/) {print $1;last;}' ?; #Collect sv_varianteffectpredictor geneBuild

    $regexp{sv_varianteffectpredictor}{assembly} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/assembly=(\S+)/) {print $1;last;}' ?; #Collect sv_varianteffectpredictor assembly

    $regexp{sv_varianteffectpredictor}{hgmd_public} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/HGMD-PUBLIC=(\S+)/) {print $1;last;}' ?; #Collect sv_varianteffectpredictor HGMD-PUBLIC version

    $regexp{sv_varianteffectpredictor}{reg_build} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/regbuild=(\S+)/) {print $1;last;}' ?; #Collect sv_varianteffectpredictor regbuild version

    $regexp{sv_varianteffectpredictor}{gencode} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/gencode=\S+\s+(\d+)/) {print $1;last;}' ?; #Collect sv_varianteffectpredictor gencode version

    $regexp{sv_vcfparser}{version} = q?perl -nae 'if($_=~/##Software=<ID=vcfParser.pl,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect sv_vcfparser version

    $regexp{sv_genmod}{version} = q?perl -nae 'if($_=~/##Software=<ID=genmod,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect SVGenmod version

    $regexp{vcftools}{version} = q?perl -nae 'if($_=~/VCFtools\s-\s(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect VCFTools version

    $regexp{plink2}{version} = q?perl -nae 'if($_=~/PLINK\s(\S+\s\S+\s\S+\s\S+\s\S+)/) {my $ret = $1;$ret =~s/\s/_/g;print $ret;last;}' ?; #Collect Plink2 version

    $regexp{variant_integrity_mendel}{fraction_of_errors}  = q?perl -nae 'unless ($_=~/^#/) {print $F[1];last;}' ?;

    $regexp{variant_integrity_mendel}{mendelian_errors}  = q?perl -nae 'unless ($_=~/^#/) {print $F[2];last;}' ?;

    $regexp{variant_integrity_father}{fraction_of_common_variants}  = q?perl -nae 'unless ($_=~/^#/) {print $F[1];last;}' ?;

    $regexp{variant_integrity_father}{common_variants}  = q?perl -nae 'unless ($_=~/^#/) {print $F[2];last;}' ?;

#$regexp{}{} = ;

    ## Writes a YAML hash to file
    write_yaml({yaml_href => \%regexp,
		yaml_file_path_ref => \$print_regexp_outfile,
	       });

}


sub help {

##help

##Function : Print help text and exit with supplied exit code
##Returns  : ""
##Arguments: $USAGE, $exit_code
##         : $USAGE    => Help text
##         : $exit_code => Exit code

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $USAGE;
    my $exit_code;

    my $tmpl = {
	USAGE => {required => 1, defined => 1, strict_type => 1, store => \$USAGE},
	exit_code => { default => 0, strict_type => 1, store => \$exit_code},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    say STDOUT $USAGE;
    exit $exit_code;
}
