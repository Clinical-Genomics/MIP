#!/usr/bin/env bash

### Creating Conda Environment: mip
conda create -y -n mip -c bioconda mosaik=2.2.26 snpeff=4.1 bcftools=1.2 samtools=1.2 bwa=0.7.12 fastqc=0.11.2 picard=1.139 bedtools=2.25 

## Activate conda environment
source activate mip 

### Install Cpanm in conda environment: mip
conda install -y -c https://conda.anaconda.org/dan_blanchard perl-app-cpanminus 

## Install Perl modules via cpanm
cpanm YAML Log::Log4perl List::MoreUtils DateTime DateTime::Format::ISO8601 DateTime::Format::HTTP DateTime::Format::Mail Set::IntervalTree LWP::Simple LWP::Protocol::https Archive::Zip DBI JSON DBD::mysql 

## Deactivate conda environment
source deactivate 

### Install PIP packages in conda environment: mip
## Activate conda environment
source activate mip 

## Install PIP packages
pip install python-Levenshtein==0.12.0 chanjo==3.0.1 cosmid==0.4.9.1 genmod==3.3.3 

## Deactivate conda environment
source deactivate 

### Install sambamba

## Create temp install directory
mkdir -p .MIP 
cd .MIP

## Download sambamba release
wget --quiet https://github.com/lomereiter/sambamba/releases/download/v0.5.9/sambamba_v0.5.9_linux.tar.bz2 -O sambamba_v0.5.9_linux.tar.bz2

## Decompress sambamba file
bzip2 -f -d sambamba_v0.5.9_linux.tar.bz2

## Extract files
tar xvf sambamba_v0.5.9_linux.tar

## Make executable
chmod 755 sambamba_v0.5.9

## Make available from conda environment
mv sambamba_v0.5.9 ~/miniconda/envs/mip/bin/

## Moving back to original working directory
cd /mnt/hds/proj/cust003/develop/modules/MIP

## Clean up
rm -rf .MIP

### Install vcfTools

## Create temp install directory
mkdir -p .MIP 
cd .MIP

## Download vcfTools
wget --quiet https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz -O vcftools-0.1.14.tar.gz

## Extract
tar xvf vcftools-0.1.14.tar.gz

## Export PERL5LIB environment variable
export PERL5LIB=/mnt/hds/proj/cust003/develop/modules/MIP/vcftools-0.1.14/src/perl/

## Move to vcfTools directory
cd vcftools-0.1.14

## Configure
./configure --prefix=/home/henrik.stranneheim/miniconda/envs/mip
make
make install

## Moving back to original working directory
cd /mnt/hds/proj/cust003/develop/modules/MIP

## Clean up
rm -rf .MIP

### Install VT

## Create temp install directory
mkdir -p .MIP 
cd .MIP

## Download VT
wget --quiet https://github.com/atks/vt/archive/0.57.tar.gz -O vt-0.57.tar.gz

## Extract
tar xvf vt-0.57.tar.gz

## Move to vt directory
cd vt-0.57

make
make test

## Make available from conda environment
mv vt ~/miniconda/envs/mip/bin/

## Moving back to original working directory
cd /mnt/hds/proj/cust003/develop/modules/MIP

## Clean up
rm -rf .MIP

### Install Plink

## Create temp install directory
mkdir -p .MIP 
cd .MIP

## Download Plink
wget --quiet http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-x86_64.zip -O plink-1.07-x86_64.zip

## Extract
unzip plink-1.07-x86_64.zip

## Move to plink directory
cd plink-1.07-x86_64

## Make available from conda environment
mv plink ~/miniconda/envs/mip/bin/

## Moving back to original working directory
cd /mnt/hds/proj/cust003/develop/modules/MIP

## Clean up
rm -rf .MIP

### Install VariantEffectPredictor

## Activate conda environment
source activate mip 

mkdir -p ~/miniconda/envs/mip/ensembl-tools-release-82/cache 

cd ~/miniconda/envs/mip

## Download VEP
wget --quiet https://github.com/Ensembl/ensembl-tools/archive/release/82.zip -O VariantEffectPredictor-82.zip

## Extract
unzip VariantEffectPredictor-82.zip

## Move to VariantEffectPredictor directory
cd ensembl-tools-release-82/scripts/variant_effect_predictor/

## Install VEP
perl INSTALL.pl --AUTO alcf -c ~/miniconda/envs/mip/ensembl-tools-release-82/cache -s homo_sapiens --ASSEMBLY GRCh37 

## Clean up
cd ~/miniconda/envs/mip

rm -rf VariantEffectPredictor-82.zip

## Deactivate conda environment
source deactivate 

