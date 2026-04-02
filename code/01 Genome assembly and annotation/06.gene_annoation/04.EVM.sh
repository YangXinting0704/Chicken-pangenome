#
#EVM
#conda install bioconda::evidencemodeler
#github:https://github.com/EVidenceModeler/EVidenceModeler

EVidenceModeler --sample_id sample --genome genome.fa --weights weights.txt --gene_predictions gene_predictions.gff3 --protein_alignments evm_pro.gff3 --segmentSize 500000 --overlapSize 50000 --CPU 120
