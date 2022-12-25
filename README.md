# Pipeline to align the reads and get genome counts from Chec-Seq experiments
The pipe.zbjline os located in **‘LAB/scripts/felix_bcl2fastq’** and consists of two scripts **dna_folder.py** and **dna_ind.py**  and a **config file** (e.g. chec_opts.txt)

**dna_folder.py:** This is the wrapping script that set the configuration for the pipeline and then submits a job for each read file in the folder to the pipeline it has the following options:

*--infolder* path to the folder containing the read files (R1 and R2)

*--outfolder* path to the output directory*--config* path to the options file (if it is inside the felix_bcl2fastq folder, the name is enough, e.g. chec_opts.txt)

*--queue wexac* queue you want to submit the individual jobs (preset: molgen-q)

*--memSize* Ram used for the pipeline (preset: 4000 MB, needs tro be increased for large files >3Mio reads)

**dna_ind.py:**This file contains the real pipeline and usually doesn’t need to be run if individual samples fail.

- First cutadapt to remove adapters

- Bowtie2 for alignmentPicard (modified) to detect pcr duplicates
- Samtools to filter and sort reads
- bedtools to align to the genome

It creates a .outfile for each readfile that contains the 5’ position of each read aligned to the genome (12 mio bases long 16Chr + ChrM)

In addition: the _stat.csv file contains the alignment and filter statistics and the config.txt contains the configuration of this run. 

**Modules necessary to run the pipeline:**
* samtools/1.8
* bedtools/2.26.0
* bowtie2/2.3.5.1
* python/3.5
* jre/8.121
* pigz/2.3.4

**Example usage:**The easiest usage of the pipeline is moving all your read files into the folder /SEQ/name_xxx and then change into this folder : cd /SEQ/name_xxx

And then just start the pipeline from inside

`python ~/../LAB/scripts/felix_bcl2fastq/dna_folder.py --memSize 8000 --queue molgen-q --infolder ./ --outfolder ./chec1`
