############################ Assembly pipeline ####################################

### In this workshop, we will reconstruct continuous (i.e. syntenic) stretches of DNA (so called contigs) from the genome of the fish Garra longipinnis (Luise's and Sandra's favourite pet's) using three types of sequencing data.

################### (1) Clone Github repository ###################

### As a first step, you will have cloned the GitHub repository of this workshop to your home directory. If this is not the case, you should do it now using the following command

git clone https://github.com/nhmvienna/Workshop_IV_DeNovoAssembly

### now, let's have a look at the data, what do we have?

cd Workshop_IV_DeNovoAssembly/data

ls -l

cd Illumina/

ls -l

### What do the raw Illumina raw data look like?

gunzip -c Garra474_1.fq.gz | head -4

## What do these top 4 rows mean? Can

### ??? Can you repeat the same for the ONT data?

################### (2) DATA Quality ###################

### Next, we will examine the data quality of the Illumina dataset using the program FASTQC

mkdir -p ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC

echo """
#!/bin/sh

## name of Job
#PBS -N fastqc

## Redirect output stream to this file.
#PBS -o ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC/fastq_log.txt

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select 10 cores and 50gb of RAM
#PBS -l select=1:ncpus=10:mem=50g

######## load dependencies #######

module load Tools/FastQC-0.11.9

## loop through all raw FASTQ and test quality

fastqc \
  --outdir ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC \
  --threads 10 \
  ~/Workshop_IV_DeNovoAssembly/data/Illumina/Garra474_1.fq.gz \
  ~/Workshop_IV_DeNovoAssembly/data/Illumina/Garra474_2.fq.gz

""" > ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC/fastqc.sh

## Submit the job to OpenPBS

qsub ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC/fastqc.sh

## check the status of your OpenPBS Job

qstat -awt

## once the job is finished, you can check the output in the browser

firefox ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC/Garra474_1_fastqc.html

## What about the ONT dataset? We will use Nanoplot for this!


mkdir -p ~/Workshop_IV_DeNovoAssembly/results/ONT_QC

echo """
#!/bin/sh

## name of Job
#PBS -N fastqc

## Redirect output stream to this file.
#PBS -o ~/Workshop_IV_DeNovoAssembly/results/ONT_QC/fastq_log.txt

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select 10 cores and 50gb of RAM
#PBS -l select=1:ncpus=10:mem=50g

######## load dependencies #######

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate nanoplot_1.32.1

######## run analyses #######

NanoPlot \
  -t 10 \
  --summary ~/Workshop_IV_DeNovoAssembly/data/ONT/sequencing_summary.txt \
  --plots dot \
  -o ~/Workshop_IV_DeNovoAssembly/results/ONT_QC

""" > ~/Workshop_IV_DeNovoAssembly/results/ONT_QC/fastqc.sh

## Submit the job to OpenPBS

qsub ~/Workshop_IV_DeNovoAssembly/results/ONT_QC/fastqc.sh

## check the status of your OpenPBS Job

qstat -awt

## once the job is finished, you can check the output in the browser

firefox ~/Workshop_IV_DeNovoAssembly/results/ONT_QC/NanoPlot-report.html

## for the PacBio data, we do not have a sequencing summary file as for the ONT data, thus, we need to use the FASTQ sequences for QC

mkdir -p ~/Workshop_IV_DeNovoAssembly/results/PacBio_QC

echo """
#!/bin/sh

## name of Job
#PBS -N fastqc

## Redirect output stream to this file.
#PBS -o ~/Workshop_IV_DeNovoAssembly/results/PacBio_QC/fastq_log.txt

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select 10 cores and 50gb of RAM
#PBS -l select=1:ncpus=10:mem=50g

######## load dependencies #######

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate nanoplot_1.32.1

######## run analyses #######

NanoPlot \
  -t 10 \
  --fastq ~/Workshop_IV_DeNovoAssembly/data/PacBio/Garra_PB.fastq.gz \
  --plots dot \
  -o ~/Workshop_IV_DeNovoAssembly/results/PacBio_QC

""" > ~/Workshop_IV_DeNovoAssembly/results/PacBio_QC/fastqc.sh

## Submit the job to OpenPBS

qsub ~/Workshop_IV_DeNovoAssembly/results/PacBio_QC/fastqc.sh

## check the status of your OpenPBS Job

qstat -awt

## once the job is finished, you can check the output in the browser

firefox ~/Workshop_IV_DeNovoAssembly/results/PacBio_QC/NanoPlot-report.html
