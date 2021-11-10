<img align="right" src="https://github.com/DavideBrex/TEs-Analyser/blob/cac80a19f5cd9fbcc70da50a9033a5b810f50279/rulegraph.svg">


# TEs-Analyser
Snakemake pipeline for Transposable Elements quantification and Differential Expression analysis


### Required packages

- DropSeq
- Bamtools
- STAR
- bedtools 
- R
- samtools 
- fastqc

### How to run

Before running, you need to set the **config.yaml** , where you need to set the paths to the STAR genome index and othe files that are required.
Also you can set there the contrast (for differential expression analysis) and the number of cpus to be used by single tools.

##### Build the container 

You can build the container available [here](https://github.com/DavideBrex/DockerFiles/tree/main/TEs-Analyser-container). Follow the istructions found in **create_cont.sh**. One you have your sif file you run the pipeline by moving to the directory whre you downloaded it and run:

        snakemake -j10 --use-singularity
        
Where 10 is the max number of cpus that snakemake will use. 
And the outputs will be stored in a directory named **results**.

##### Author
Davide Bressan - Fulvio Chiacchiera Lab

##### Disclaimer
Some parts of the code where taken from the amazing work of Daniel Fernandez Perez. Specifically, from his repository: [Rna-seq pipeline](https://github.com/dfernandezperez/RNAseq-Snakemake)
