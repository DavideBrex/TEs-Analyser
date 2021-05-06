import pandas as pd
import numpy as np
import yaml
from snakemake.utils import validate, min_version


##### set minimum snakemake version #####
min_version("5.4.3")

##### Singularity image path
#singularity: "/hpcnfs/data/DP/Singularity/dfernandezperez-bioinformatics-singularity-master-chipseq.simg"

##### config file
configfile: "configuration/config.yaml"

#read samples info
SAMPLES=pd.read_csv(config['samples'], sep = "\t").set_index("NAME", drop=False).sort_index()
#read fastq files
FASTQ_FILES=pd.read_csv(config["fastq_files"], dtype=str, sep ="\t").set_index(["sample", "lane"], drop=False).sort_index()


print(SAMPLE)