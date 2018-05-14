# Setup

## Fastq filename convention
The permanent filename should follow the following format:

``{LANE}_{DATE}_{FLOW-CELL}_{SAMPLE-ID}_{BARCODE-SEQ}_{DIRECTION 1/2}.fastq[.qz]``

Where some types or formats are required some each element:
- LANE = Integer
- DATE = YYMMDD
- BARCODE-SEQ = A, C, G, T or integer
- DIRECTION = 1 or 2

The `familyID` and `sampleID(s)` needs to be unique and the sample id supplied should be equal to the {SAMPLE_ID} in the filename.
Underscore cannot be part of any element in the file name as this is used as the seperator for each element.

However, MIP will except filenames in other formats as long as the filename contains the sample id and the mandatory information can be collected from the fastq header.

## Meta-Data
MIP requires pedigree information recorded in a pedigree.yaml file and a config file.

* [Pedigree file] \(YAML-format\)
* [Configuration file] \(YAML-format\)

## Dependencies
Make sure you have installed all dependencies and that they are in your ``$PATH``.
You only need to install the dependencies that are required for the modules that you want to run. If you have not installed a dependency for a module, MIP will tell you what dependencies you need to install (or add to your ``$PATH``) and exit. MIP comes with an install script, which will install all necessary programs to execute models in MIP via bioconda and/or $SHELL.

### **Programs**

- Simple Linux Utility for Resource Management ([SLURM]) (version: 2.6.0)
- [Bcftools] (version: 1.6)
- [BedTools] (version: 2.26.0)
- [BWA] (version: 0.7.15)
- [BWAKit] (version: 0.7.12)
- [Chanjo] (version: 4.2.0)
- [Cnvnator] (version: 0.3.3)
- [Delly] (version: 0.7.7)
- [FastQC] (version: 0.11.5)
- [Freebayes] (version: 1.1.0)
- [GATK] (version: 3.8)
- [GENMOD] (version: 3.7.2)
- [Htslib] (version: 1.6)
- [Manta] (version: 1.1.0)
- [MultiQC] (version: 1.4)
- [Peddy] (version: 0.3.1)
- [PicardTools] (version: 2.14.1)
- [PLINK2] (version: 1.90b3x35)
- [Sambamba] (version: 0.6.6)
- [Samtools] (version: 1.6)
- [SnpEff] (version: 4.3.1)
- [Svdb] (version: 1.0.7)
- [Tiddit] (version: 2.2.1)
- [Variant_integrity] (version: 0.0.4)
- [Vcf2cytosure] (version: 0.2.0)
- [Vcfanno] (version: 0.1.0)
- [VEP] (version: 91) with plugin "MaxEntScan, LoFtool"
- [VT] (version: 20151110)

The version number after the software name are tested for compatibility with MIP.

### Databases/References

MIP can build/download many program prerequisites automatically via the mip_install script using flag ``--reference_dir [reference_dir]``, which will use the MIP script ``download_reference.pl``.

### **Automatic Build:**

Human Genome Reference Meta Files:
 1. The sequence dictionnary (".dict")
 2. The ".fasta.fai" file

BWA:
 1. The BWA index of the human genome.

#### *Note*
If you do not supply these parameters (Bwa) MIP will create these from scratch using the supplied human reference genom as template.

Capture target files:
 1. The "infile_list" and .pad100.infile_list files used in ``ppicardtoolscollecthsmetrics``.
 2. The ".pad100.interval_list" file used by some GATK modules.

#### *Note*
If you do not supply these parameters MIP will create these from scratch using the supplied "latest" supported capture kit ".bed" file and the supplied human reference genome as template.

[Bcftools]: http://www.htslib.org/
[BedTools]: http://bedtools.readthedocs.org/en/latest/
[BWA]: https://github.com/lh3/bwa
[BWAKit]: https://github.com/lh3/bwa/tree/master/bwakit
[Chanjo]: https://chanjo.readthedocs.org/en/latest/
[Cnvnator]: https://github.com/abyzovlab/CNVnator
[Configuration file]: https://github.com/henrikstranneheim/MIP/blob/master/templates/mip_config.yaml
[Delly]: https://github.com/dellytools/delly/
[FastQC]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[Freebayes]: https://github.com/ekg/freebayes
[GATK]: http://www.broadinstitute.org/gatk/
[GENMOD]: https://github.com/moonso/genmod/
[Htslib]: http://www.htslib.org/
[Manta]: https://github.com/Illumina/manta
[MultiQC]: https://github.com/ewels/MultiQC
[Peddy]: https://github.com/brentp/peddy
[Pedigree file]: https://github.com/Clinical-Genomics/MIP/tree/master/templates/643594-miptest_pedigree.yaml   
[PicardTools]: http://broadinstitute.github.io/picard/
[PLINK2]: https://www.cog-genomics.org/plink2
[Sambamba]: http://lomereiter.github.io/sambamba/
[Samtools]: http://www.htslib.org/
[SLURM]: http://slurm.schedmd.com/
[SnpEff]: http://snpeff.sourceforge.net/
[Svdb]: https://github.com/J35P312/SVDB
[Tabix]: http://samtools.sourceforge.net/tabix.shtml
[Tiddit]: https://github.com/J35P312/TIDDIT
[Variant_integrity]: https://github.com/moonso/variant_integrity
[Vcf2cytosure]: https://github.com/NBISweden/vcf2cytosure
[Vcfanno]: https://github.com/brentp/vcfanno
[VEP]: https://github.com/Ensembl/ensembl-vep
[VT]: https://github.com/atks/vt
