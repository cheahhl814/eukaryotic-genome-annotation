# BUSCO — busco.ezlab.org (fetched 2026-10-21)

BUSCO metric is complementary to technical metrics like N50; based on near-universal single-copy orthologs (OrthoDBv12.2 datasets).

## Obtain BUSCO
- Available as a conda package (bioconda `busco`) and as a Docker image; installation guide at busco.ezlab.org/userguide.
- Display all available datasets: `busco --list-datasets`; manual download: https://busco-data.ezlab.org/v6/data/lineages/

## Mandatory options
```
-i  # Input sequence file or folder
-m  # BUSCO analysis mode: 'genome', 'protein' or 'transcriptome'
```

## Recommended options
```
-c  # Number of threads/cores
-l  # BUSCO lineage dataset for scoring
-o  # Output folder name
```

## Pipeline-specific options
```
--augustus             # BLAST/Augustus pipeline
--augustus_species     # species for Augustus training
--augustus_parameters  # "--PARAM1=VALUE1,--PARAM2=VALUE2"
--metaeuk              # Metaeuk pipeline
--metaeuk_parameters, --metaeuk_rerun_parameters
--miniprot             # miniprot pipeline (default for eukaryota)
```

## Other options
```
--config CONFIG_FILE   # config file
--download [dataset]   # dataset name, all, prokaryota, eukaryota, virus
--download_path PATH
-f  # force overwrite output dir
-r  # restart incomplete run
--list-datasets
--offline  # no attempt to download
-q  # quiet
--tar  # compress output subdirs
--skip_bbtools
--scaffold_composition
--contig_break n       # default n=10
```

## Citations
- odb12 datasets: doi 10.1093/nar/gkae987; odb10: doi 10.1093/molbev/msab199
