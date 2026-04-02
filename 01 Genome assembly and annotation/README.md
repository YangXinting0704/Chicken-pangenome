Note: some codes were adopted according to the configuration environment of the Institute of Animal Sciences,Chinese Academy of Agricultural Sciences.

## Description for each directory:

- The [01.draft assembly] directory contains codes for genome assembly and genome assembly flowchart. 

- The [02.chromosome anchoring] directory contains codes for eliminate redundant sequences，shord reads mapping, hic reads mapping, chromosome anchoring.

- The [03.gap filling and telomere filling] directory contains codes for genome alignment,spanning the gap regions,gap filling and resolve missing telomeres.

- The [04.assessment of assembly quality] directory contains codes for Statistical assembly results,BUSCO evaluation,Measure the completeness of the genome and the quality of the assembly,etc.The config file for HiC-Pro is also provided.

- The [05.repeat annotation] directory contains codes for TE annotation,Construction of de novo repeat library.

- The [06.gene_annoation] directory contains codes for protein-coding gene structure,ab initio prediction, transcriptome-based prediction, and homology-based prediction, integrate all gene model evidences,functional annotations.The weight file for EVM also is also provided.

- The [07.centromere identification] direcotory contains codes for  chip-seq analysis for CENP-A, visualization tandem repeats using [StainedGlass].The config.yaml for StainedGlass also is also provided.

- The [08.methylation determination] directory contains codes for Methylation determination using ONT data.






