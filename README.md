<img align="right" src="https://raw.githubusercontent.com/DavideBrex/TEs-Analyser/67e1b375aebb02d3a3cc4a8da651d7a698b85c3c/rulegraph_2.svg">


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


You can build the container available [here](https://github.com/DavideBrex/DockerFiles/tree/main/TEs-Analyser-container). Follow the istructions found in **create_cont.sh**. One you have your sif file you run the pipeline by moving to the directory whre you downloaded it and run:

        snakemake -j10 --use-singularity
        
And the outputs will be stored in a directory named **results**.

##### Author
Davide Bressan - Fulvio Chiacchiera Lab

##### Disclaimer
Some parts of the code where taken from the amazing work of Daniel Fernandez Perez. Specifically, from his repository: [Rna-seq pipeline](https://github.com/dfernandezperez/RNAseq-Snakemake)
