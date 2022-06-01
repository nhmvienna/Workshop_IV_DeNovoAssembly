#### Assembly pipeline ####

### In this workshop, we will to reconstruct continuous (i.e. syntenic) stretches of DNA (so called contigs) from the genome of the fish Garra longipinnis (Luise's and Sandra's favourite pet's) using three types of sequencing data.

### (1) Clone Github repository

### As a first step, you will have cloned the GitHub repository of this workshop to your home directory. If this is not the case, you should do it now using the following command

git clone https://github.com/nhmvienna/Workshop_IV_DeNovoAssembly

### now, let's have a look at the data, what do we have?

cd Workshop_IV_DeNovoAssembly/data/

ls -l

cd Illumina/

ls -l

### What do the raw Illumina raw data look like?

gunzip -c Garra474_1.fq.gz | head -4

## What do these top 4 rows
