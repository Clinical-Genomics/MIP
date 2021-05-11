# Setup

## Fastq filename convention
The permanent filename should follow the following format:

``{LANE}_{DATE}_{FLOW-CELL}_{SAMPLE-ID}_{BARCODE-SEQ}_{DIRECTION 1/2}.fastq[.qz]``

Where some types or formats are required for each element:
- LANE = Integer
- DATE = YYMMDD
- BARCODE-SEQ = A, C, G, T or integer
- DIRECTION = 1 or 2

The `case_id` and `sample_id(s)` needs to be unique and the sample id supplied should be equal to the {SAMPLE_ID} in the filename.
Underscore cannot be part of any element in the file name as this is used as the separator for each element.

However, MIP will accept filenames in other formats as long as the filename contains the sample id and the mandatory information can be collected from the fastq header within the file.

## Meta-Data
MIP requires pedigree information recorded in a pedigree.yaml file and a config file.

* [Pedigree file] \(YAML-format\)
* [Configuration file] \(YAML-format\)

## Dependencies
MIP comes with an install application, which will install all necessary programs to execute models in MIP via conda and/or $SHELL. Make sure you have installed all dependencies via the MIP install application and that you have loaded your MIP base environment.
You only need to install the dependencies that are required for the recipes that you want to run. If you have not installed a dependency for a module, MIP will tell you what dependencies you need to install and exit.

**Extra CPANM modules**
You can speed up, for instance, the Readonly module by also installing the companion module Readonly::XS. No change to the code is required and the Readonly module will call the Readonly::XS module if available.  

### **Programs**

- Simple Linux Utility for Resource Management ([SLURM]) (version: 18.08.0)

#### **Pipeline: Rare disease**
- [Arriba] (version: 2.1.0)
- [Bcftools] (version: 1.11)
- [BWA] (version: 0.7.17)
- [BedTools] (version: 2.30.0)
- [Blobfish] (version: 0.0.2 )
- [BootstrapAnn] (version: df02f35t)
- [BWAKit] (version: 0.7.17)
- [CADD] (version: 1.6)
- [Chanjo] (version: 4.6.0)
- [Chromograph] (version: 1.1)
- [Cnvnator] (version: 0.4.1)
- [Cyrius] (version: v1.1)
- [Expansionhunter] (version 4.0.2)
- [FastQC] (version: 0.11.9)
- [Deepvariant] (version: 1.1.0)
- [Delly] (version 0.8.7)
- [GATK] (version: 3.8.1 and 4.2.0.0)
- [GENMOD] (version: 3.7.3)
- [Gffcompare] (version: 0.11.2)
- [Glnexus] (version: 1.3.1)
- [Htslib] (version: 1.12)
- [Manta] (version: 1.6.0)
- [Megafusion] (version: cc6676c)
- [MultiQC] (version: 1.10.1)
- [Pdfmerger] (version: 1.0)
- [Peddy] (version: 0.4.3)
- [PicardTools] (version: 2.25.0)
- [PLINK] (version: 1.90b3x35)
- [Preseq] (version: 3.1.2)
- [Rhocall] (version: 0.5.1)
- [RSeQC] (version: 4.0.0)
- [rtg-tools] (version: 3.12)
- [Salmon] (version: 1.4.0)
- [Sambamba] (version: 0.6.8)
- [Samtools] (version: 1.11)
- [SMNCopyNumberCaller] (version: v1.1.1)
- [STAR-Fusion] (version: 1.10.0)
- [STAR] (version: 2.7.8a)
- [Stranger] (version: 0.7)
- [StringTie] (version: 2.1.3b)
- [Svdb] (version: 2.4.0)
- [Telomerecat] (version: 3.4.0)
- [Tiddit] (version: 2.12.1)
- [Upd] (version: 0.1.1)
- [Varg] (version: 1.2.0)
- [Vcf2cytosure] (version: 0.5.1)
- [Vcfanno] (version: 0.3.2)
- [VEP] (version: 103.1) with plugin "ExACpLI", "MaxEntScan, LoFtool", "SpliceAI"

The version number after the software name are tested for compatibility with MIP.

### Databases/References

MIP can download many program prerequisites automatically via the mip download application ``mip download [PIPELINE]``.

MIP will build references and meta files (if required) prior to starting an analysis pipeline ``mip analyse [PIPELINE]``.

### **Automatic Build:**

Human Genome Reference Meta Files:
 1. The sequence dictionary (".dict")
 2. The ".fasta.fai" file

BWA:
 1. The BWA index of the human genome.

Star:
 1. Star index files of the human genome

#### *Note*
If you do not supply these parameters (Bwa/Star) MIP will create these from scratch using the supplied human reference genome as template.

Capture target files:
 1. The "infile_list" and .pad100.infile_list files used in ``picardtools_collecthsmetrics``.
 2. The ".pad100.interval_list" file used by some GATK recipes.

#### *Note*
If you do not supply these parameters MIP will create these from scratch using the supplied "latest" supported capture kit ".bed" file and the supplied human reference genome as template.

### Private References
Some references are not available for download because they require a license or contain data that are not consented for sharing. These references have to be manually applied for and added to the analysis where appropriate:

#### SweFreq — The Swedish Frequency resource for genomics
This dataset contains whole-genome variant frequencies for 1000 Swedish individuals generated within the SweGen project. One can request data access and download files from: https://swefreq.nbis.se/

Corresponding MIP references:
 - grch37_anon-swegen_str_nsphs_-1000samples-.vcf.gz (Autozygosity calculation;Rhocall)
 - grch37_anon_swegen_snp_-2016-10-19-.tab.gz (Frequency annotation;Vcfanno)
 - grch37_anon-swegen_indel_-1000samples-.vcf.gz (Frequency annotation;Vcfanno)
 - grch37_swegen_concat_sort_-20170830-.vcf (Structural variant frequency annotation; Svdb)

#### Spidex - Splicing prediction
SPIDEX is a computational model that uses the Percentage of Spliced-In (PSI) metric to evaluate whether a certain splicing isoform is more enriched under the presence/absence of a given variant. Unfortunately, the Deepgenomics company that used to provide spidex scores seem to have shut down. The files can be downloaded via Annovar, which however also requires a license.

Corresponding MIP references:
 - grch37_spidex_public_noncommercial_-v1_0-.tab.gz (Splicing annotation; Genmod)

#### Local frequency Databases
We use several local frequency databases, that we unfortunately are not allowed to share, but can be built using locusdb or Svdb: https://github.com/moonso/loqusdb.

Corresponding MIP references:
 - grch37_loqusdb_snv_indel_-2018-12-18-.vcf.gz (SNV/INDELS; Vcfanno)
 - grch37_mip_sv_svdb_export_-2018-10-09-.vcf (SV; Svdb)
 - grch37_svdb_query_clingen_ngi_-v1.0.0-.vcf (SV;Svdb)
 - grch37_svdb_query_decipher_-v1.0.0-.vcf (local array SVs frequency annotation; Svdb)

#### Local clincial significance databases
 Variants annotated as benign or pathogenic from array data.

 Corresponding MIP references:
 - grch37_svdb_query_clingen_cgh_benign_-v1.0.0-.vcf
 - grch37_svdb_query_clingen_cgh_pathogenic_-v1.0.0-.vcf

#### GATK exome dataset
GATK needs more data in the variant calling for exomes than a single sample or trio. MIP adds in other previously sequenced samples in the variant calling as a supplementary dataset.

Corresponding MIP references:
 - grch37_gatk_merged_reference_samples.txt

## Gene panel for the clinical test
MIP will split the variants into two sets (clinical a.k.a "selected" and research) based on gene coordinates and hgnc_id, which is recorded in gene_panels.bed file(s) using MIPs own vcfparser. A template for grch37 can be found in MIPs dir under `templates/gene_panels.bed`.  

#### SpliceAI - Splicing predicition
Annotations for all possible substitutions, 1 base insertions and 1-4 base deletions. Available from https://basespace.illumina.com/s/5u6ThOblecrh

Corresponding MIP references:
 - grch37_spliceai_scores_raw_snv_-v1.3-.vcf.gz
 - grch37_spliceai_scores_raw_indel_-v1.3-.vcf.gz

[Arriba]: https://github.com/suhrig/arriba
[Bcftools]: http://www.htslib.org/
[BedTools]: http://bedtools.readthedocs.org/en/latest/
[Blobfish]: https://github.com/J35P312/BlobFish
[BoostrapAnn]: https://github.com/J35P312/BootstrapAnn
[BWA]: https://github.com/lh3/bwa
[BWAKit]: https://github.com/lh3/bwa/tree/master/bwakit
[CADD]: (https://github.com/kircherlab/CADD-scripts)
[Chanjo]: https://chanjo.readthedocs.org/en/latest/
[Chromograph]: https://github.com/mikaell/chromograph
[Cnvnator]: https://github.com/abyzovlab/CNVnator
[Configuration file]: ../templates/mip_rd_dna_config.yaml
[Cyrius]: https://github.com/Illumina/Cyrius
[Delly]: https://github.com/dellytools/delly
[Deepvariant]: https://github.com/google/deepvariant
[Expansionhunter]: https://github.com/Illumina/ExpansionHunter
[FastQC]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[GATK]: http://www.broadinstitute.org/gatk/
[GENMOD]: https://github.com/moonso/genmod/
[Gffcompare]: https://ccb.jhu.edu/software/stringtie/gffcompare.shtml
[Glnexus]: https://github.com/dnanexus-rnd/GLnexus
[Htslib]: http://www.htslib.org/
[Manta]: https://github.com/Illumina/manta
[Megafusion]: https://github.com/J35P312/MegaFusion
[MultiQC]: https://github.com/ewels/MultiQC
[Pdfmerger]: https://github.com/northwestwitch/pdfmerger
[Peddy]: https://github.com/brentp/peddy
[Pedigree file]: ../templates/643594-miptest_pedigree.yaml  
[PicardTools]: http://broadinstitute.github.io/picard/
[PLINK]: https://www.cog-genomics.org/plink2
[Preseq]: http://smithlabresearch.org/software/preseq/
[Rhocall]: https://github.com/dnil/rhocall
[RSeQC]: http://rseqc.sourceforge.net/
[rtg-tools]: https://github.com/RealTimeGenomics/rtg-tools
[Salmon]: https://combine-lab.github.io/salmon/
[Sambamba]: http://lomereiter.github.io/sambamba/
[Samtools]: http://www.htslib.org/
[SLURM]: http://slurm.schedmd.com/
[SMNCopyNumberCaller]: https://github.com/Illumina/SMNCopyNumberCaller
[STAR]: https://github.com/alexdobin/STAR
[STAR-Fusion]: https://github.com/STAR-Fusion/STAR-Fusion
[Stranger]: https://github.com/moonso/stranger
[StringTie]: https://ccb.jhu.edu/software/stringtie/
[Svdb]: https://github.com/J35P312/SVDB
[Tabix]: http://samtools.sourceforge.net/tabix.shtml
[Telomerecat]: https://github.com/cancerit/telomerecat
[Tiddit]: https://github.com/J35P312/TIDDIT
[Upd]: https://github.com/bjhall/upd
[Varg]: https://github.com/Clinical-Genomics/varg  
[Vcf2cytosure]: https://github.com/NBISweden/vcf2cytosure
[Vcfanno]: https://github.com/brentp/vcfanno
[VEP]: https://github.com/Ensembl/ensembl-vep
