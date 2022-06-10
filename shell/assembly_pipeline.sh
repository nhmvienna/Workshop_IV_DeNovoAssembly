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

gunzip -c Garra_Ill_R1.fq.gz | head -4

## What do these top 4 rows mean?

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
  ~/Workshop_IV_DeNovoAssembly/data/Illumina/Garra_Ill_R1.fq.gz \
  ~/Workshop_IV_DeNovoAssembly/data/Illumina/Garra_Ill_R2.fq.gz

""" > ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC/fastqc.sh

## Submit the job to OpenPBS

qsub ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC/fastqc.sh

## check the status of your OpenPBS Job

qstat -awt

## once the job is finished, you can check the output in the browser

firefox ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC/Garra_Ill_R1_fastqc.html
firefox ~/Workshop_IV_DeNovoAssembly/results/Illumina_QC/Garra_Ill_R2_fastqc.html

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

""" > ~/Workshop_IV_DeNovoAssembly/results/ONT_QC/nanoplot_ont.sh

## Submit the job to OpenPBS

qsub ~/Workshop_IV_DeNovoAssembly/results/ONT_QC/nanoplot_ont.sh

## check the status of your OpenPBS Job

qstat -awt

## once the job is finished, you can check the output in the browser

firefox ~/Workshop_IV_DeNovoAssembly/results/ONT_QC/NanoPlot-report.html

################### (3) Trimming of Illumina reads ###################

## Before we start the actual assembly, we need to "clean up" the Illumina reads, i.e. to trim away tails of the reads with low quality and adaptor sequences that were used for Illumina sequencing. We use the program trim_galore for that.

mkdir ~/Workshop_IV_DeNovoAssembly/results/trimmed

echo """
#!/bin/sh

## name of Job
#PBS -N trim_galore

## Redirect output stream to this file.
#PBS -o ~/Workshop_IV_DeNovoAssembly/results/trimmed/log.txt

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select 10 cores and 50gb of RAM
#PBS -l select=1:ncpus=10:mem=50g

######## load dependencies #######

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate trim-galore-0.6.2

## loop through all FASTQ pairs and trim by quality PHRED 20, min length 85bp and automatically detect & remove adapters

cd ~/Workshop_IV_DeNovoAssembly/results/trimmed

trim_galore \
  --paired \
  --quality 20 \
  --length 85  \
  --cores 10 \
  --fastqc \
  --gzip \
  ~/Workshop_IV_DeNovoAssembly/data/Illumina/Garra_Ill_R1.fq.gz \
  ~/Workshop_IV_DeNovoAssembly/data/Illumina/Garra_Ill_R2.fq.gz

""" > ~/Workshop_IV_DeNovoAssembly/results/trimmed/trim.sh

qsub ~/Workshop_IV_DeNovoAssembly/results/trimmed/trim.sh

## check the status of your OpenPBS Job
qstat -awt

## once the job is finished, you can check the quality of the trimmed reads in the browser
firefox ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R1_val_1_fastqc.html
firefox ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R2_val_2_fastqc.html

################### (4) Genome-size estimation ###################

## Now that we have an idea the quality and trimmed Illumina reads we can use these data to get a rough idea about the expected size of the genome. Note that we are using only a very small subset. Thus, the estimate will only be very rough.

## First, we will use the program JellyFish to count the number of k-mers in the dataset. A k-mer is a unqiue sequence of a given length n (for example n=31bp) found in the pool of sequences. The original reads will therefore chopped down into substrings of size n, for example (n=5):

# ACGGTGAGGAT
# ACGGT
#  CGGTG
#   GGTGA
#    GTGAG
#     TGAGG
#      GAGGA
#       AGGAT

## After that, the Program GenomeScope will estimate the Genome-size based on the coverage distribution of the k-mers. See https://github.com/nhmvienna/AutDeNovo#3-genome-size-estimation for more details.

mkdir ~/Workshop_IV_DeNovoAssembly/results/genomesize

echo """
  #!/bin/sh

  ## name of Job
  #PBS -N jellyfish

  ## Redirect output stream to this file.
  #PBS -o ~/Workshop_IV_DeNovoAssembly/results/genomesize/jellyfish_log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select a maximum walltime of 2h
  #PBS -l walltime=48:00:00

  ## Select 10 cores and 50gb of RAM
  #PBS -l select=1:ncpus=10:mem=50gb

  ## load all necessary software into environment
  module load Assembly/Jellyfish-2.3.0
  module load Assembly/genomescope-2.0

  ## unzip files
  gunzip -c ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R1_val_1.fq.gz \
  > ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R1_val_1.fq &
  gunzip -c ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R2_val_2.fq.gz \
  > ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R2_val_2.fq

  wait

  ## run Jellyfish
  ## parameters
  # -C canononical; count both strands
  # -m 31 Length of mers
  # -s initial hash size

  jellyfish-linux count \
    -C \
    -m 31 \
    -s 100M \
    -t 10 \
    -F 2 \
    -o ~/Workshop_IV_DeNovoAssembly/results/genomesize/reads.jf \
    ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R1_val_1.fq \
    ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R2_val_2.fq

  ## remove unzipped copy of reads
  rm -f ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R*_val_*.fq

  ## make a histogram of all k-mers
  jellyfish-linux histo \
    -t 10 \
    ~/Workshop_IV_DeNovoAssembly/results/genomesize/reads.jf \
    > ~/Workshop_IV_DeNovoAssembly/results/genomesize/reads.histo

  ## run GenomeScope
  genomescope.R \
  -i ~/Workshop_IV_DeNovoAssembly/results/genomesize/reads.histo \
  -k 31 \
  -p 2 \
  -o ~/Workshop_IV_DeNovoAssembly/results/genomesize/stats

""" > ~/Workshop_IV_DeNovoAssembly/results/genomesize/genomesize.sh

qsub  ~/Workshop_IV_DeNovoAssembly/results/genomesize/genomesize.sh

################### (5) De Novo Assembly ###################

## Now that we have an idea about the Approximate genome-size we can start the de-novo assembly

## When using Illumina data alone or in combination with single molecule sequencing data (ONT and PacBio), we will use SPAdes which reconstructs genomes using DeBrujin-Graph algorithms. If only single molecule sequencing data is available, we will use the FLYE assembler, which is based on Repeat Graphs

## First, we will use SPAdes with Illumina and ONT data. In this case, the ONT data will be used to combine contigs (from Illumina data) into scaffolds.

mkdir -p ~/Workshop_IV_DeNovoAssembly/results/denovo/spades

echo """
  #!/bin/sh

  ## name of Job
  #PBS -N denovo_Spades

  ## Redirect output stream to this file.
  #PBS -o ~/Workshop_IV_DeNovoAssembly/results/denovo/spades/log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select a maximum walltime of 2h
  #PBS -l walltime=100:00:00

  ## Select 10 cores and 50gb of RAM
  #PBS -l select=1:ncpus=10:mem=50gb

  ## load all necessary software into environment
  module load Assembly/SPAdes_3.15.3

  ## first concatenate all ONT reads into a single file
  cat ~/Workshop_IV_DeNovoAssembly/data/ONT/Garra_ONT_*.fastq.gz \
  > ~/Workshop_IV_DeNovoAssembly/data/ONT/Garra_ONT.fastq.gz

  spades.py \
    -1 ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R1_val_1.fq.gz \
    -2 ~/Workshop_IV_DeNovoAssembly/results/trimmed/Garra_Ill_R2_val_2.fq.gz \
    --nanopore ~/Workshop_IV_DeNovoAssembly/data/ONT/Garra_ONT.fastq.gz \
    -t 10 \
    -m 50 \
    -o ~/Workshop_IV_DeNovoAssembly/results/denovo/spades

  rm -f ~/Workshop_IV_DeNovoAssembly/data/ONT/Garra_ONT.fastq.gz

""" > ~/Workshop_IV_DeNovoAssembly/results/denovo/spades/spades.sh

qsub ~/Workshop_IV_DeNovoAssembly/results/denovo/spades/spades.sh

## alternatively, we also try FLYE with the ONT data only

mkdir -p ~/Workshop_IV_DeNovoAssembly/results/denovo/flye

echo """
  #!/bin/sh

  ## name of Job
  #PBS -N denovo_Flye

  ## Redirect output stream to this file.
  #PBS -o ~/Workshop_IV_DeNovoAssembly/results/denovo/flye/log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select a maximum walltime of 2h
  #PBS -l walltime=100:00:00

  ## Select 10 cores and 50gb of RAM
  #PBS -l select=1:ncpus=10:mem=50gb

  ## load all necessary software into environment
  source /opt/anaconda3/etc/profile.d/conda.sh
  conda activate flye-2.9

  ## first concatenate all ONT reads into a single file
  cat ~/Workshop_IV_DeNovoAssembly/data/ONT/Garra_ONT_*.fastq.gz \
  > ~/Workshop_IV_DeNovoAssembly/data/ONT/Garra_ONT.fastq.gz

  flye \
  --nano-raw ~/Workshop_IV_DeNovoAssembly/data/ONT/Garra_ONT.fastq.gz \
  --out-dir ~/Workshop_IV_DeNovoAssembly/results/denovo/flye \
  --threads 10 \
  --scaffold

""" > ~/Workshop_IV_DeNovoAssembly/results/denovo/flye/flye.sh

qsub ~/Workshop_IV_DeNovoAssembly/results/denovo/flye/flye.sh

#### OK, now the assemblies is finished, what now?

################### (6) Assembly statistics ###################

## At first, we will produce some statistics that summarize the assembly. These include, for example, (1) the total number of contigs, (2) the N50 value, i.e. N50 is defined as the sequence length of the shortest contig at 50% of the total genome length (after sorting the contigs by decreasing length). We use the program QUAST to calculate these

mkdir -p ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades/

echo """

  #!/bin/sh

  ## name of Job
  #PBS -N QUAST_Spades

  ## Redirect output stream to this file.
  #PBS -o~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades/log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select 10 cores and 50gb of RAM
  #PBS -l select=1:ncpus=10:mem=50g

  ######## load dependencies #######

  module load Assembly/Quast-5.1.0rc1

  ######## run analyses #######

  quast.py \
  --output-dir ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades \
  --threads 10 \
  --eukaryote \
  -f \
  ~/Workshop_IV_DeNovoAssembly/results/denovo/spades/scaffolds.fasta

""" > ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades/quast.sh

qsub ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades/quast.sh

### now we repeat the same for the FLYE assembly

mkdir -p ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye/

echo """

  #!/bin/sh

  ## name of Job
  #PBS -N QUAST_Flye

  ## Redirect output stream to this file.
  #PBS -o~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye/log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select 10 cores and 50gb of RAM
  #PBS -l select=1:ncpus=10:mem=50g

  ######## load dependencies #######

  module load Assembly/Quast-5.1.0rc1

  ######## run analyses #######

  quast.py \
  --output-dir ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye \
  --threads 10 \
  --eukaryote \
  -f \
  ~/Workshop_IV_DeNovoAssembly/results/denovo/flye/assembly.fasta

""" > ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye/quast.sh

qsub ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye/quast.sh

### last, but not least, we will test how many BUSCO (XXX) genes can be detected in our de novo assembly

## First, for the Illumina-based assembly

mkdir ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades/Busco

echo """

  #!/bin/sh

  ## name of Job
  #PBS -N BUSCO_flye

  ## Redirect output stream to this file.
  #PBS -o ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades/Busco/log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select 10 cores and 50gb of RAM
  #PBS -l select=1:ncpus=10:mem=50g

  ######## load dependencies #######

  source /opt/anaconda3/etc/profile.d/conda.sh
  conda activate busco_5.2.2

  ######## run analyses #######

  ## Go to pwd
  cd ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades/Busco

  busco -i ~/Workshop_IV_DeNovoAssembly/results/denovo/spades/scaffolds.fasta \
    -o spades \
    -m genome \
    -c 10 \
    -f \
    -l vertebrata_odb10

""" > ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades/Busco/Spades_busco.sh

qsub ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/spades/Busco/Spades_busco.sh

## ...and then for the ONT-based assembly

mkdir ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye/Busco

echo """

  #!/bin/sh

  ## name of Job
  #PBS -N BUSCO_flye

  ## Redirect output stream to this file.
  #PBS -o ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye/Busco/log.txt

  ## Stream Standard Output AND Standard Error to outputfile (see above)
  #PBS -j oe

  ## Select 10 cores and 50gb of RAM
  #PBS -l select=1:ncpus=10:mem=50g

  ######## load dependencies #######

  source /opt/anaconda3/etc/profile.d/conda.sh
  conda activate busco_5.2.2

  ######## run analyses #######

  ## Go to pwd
  cd ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye/Busco

  busco -i ~/Workshop_IV_DeNovoAssembly/results/denovo/flye/assembly.fasta \
    -o Flye \
    -m genome \
    -c 10 \
    -f \
    -l vertebrata_odb10

""" > ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye/Busco/Flye.sh

qsub ~/Workshop_IV_DeNovoAssembly/results/AssemblyQC/flye/Busco/Flye.sh

### OK; lot's of work! wouldn't it be nice to do everything in one go?

mkdir ~/Workshop_IV_DeNovoAssembly/results/Automated_ILL

sh /media/inter/pipelines/AutDeNovo/AutDeNovo.sh

## Let's try it out with our test Illumina dataset

sh /media/inter/pipelines/AutDeNovo/AutDeNovo.sh \
  Name=ILL \
  OutputFolder=~/Workshop_IV_DeNovoAssembly/results/Automated_ILL \
  Fwd=~/Workshop_IV_DeNovoAssembly/data/Illumina/Garra_Ill_R1.fq.gz \
  Rev=~/Workshop_IV_DeNovoAssembly/data/Illumina/Garra_Ill_R1.fq.gz \
  threads=10 \
  RAM=20 \
  RAMAssembly=20 \
  decont=no \
  SmudgePlot=no \
  BuscoDB=vertebrata_odb10

## Now repeat with the ONT data

mkdir ~/Workshop_IV_DeNovoAssembly/results/Automated_ONT

sh /media/inter/pipelines/AutDeNovo/AutDeNovo.sh \
  Name=ONT \
  OutputFolder=~/Workshop_IV_DeNovoAssembly/results/Automated_ONT \
  ONT=~/Workshop_IV_DeNovoAssembly/data/ONT \
  threads=10 \
  RAM=20 \
  RAMAssembly=20 \
  decont=no \
  SmudgePlot=no \
  BuscoDB=vertebrata_odb10

## what do you get??
