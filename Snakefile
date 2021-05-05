import pandas as pd
import numpy as np
import yaml
from snakemake.utils import validate, min_version


##### set minimum snakemake version #####
min_version("5.4.3")

##### Singularity image path
#singularity: "/hpcnfs/data/DP/Singularity/dfernandezperez-bioinformatics-singularity-master-chipseq.simg"

##### config file
#configfile: "configuration/config.yaml"