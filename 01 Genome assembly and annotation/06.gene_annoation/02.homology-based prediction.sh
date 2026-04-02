#homology-based prediction
#!/bin/bash
#est_isoseq
#conda install bioconda::transdecoder
#conda install bioconda::stringtie
#conda install bioconda::blast
#conda install bioconda::hmmer
#github:https://github.com/TransDecoder/TransDecoder https://github.com/ncbi/blast_plus_docs https://github.com/EddyRivasLab/hmmer

minimap2 -I 30g  -ax splice -f 1000 --sam-hit-only --secondary=no --eqx -K 100M -t 40 --cap-sw-mem=3g genome.fa  est_isoseq.fa | samtools view -Sb - > est_isoseq.bam
samtools sort -@ 20 -o sort.bam est_isoseq.bam
samtools index  sort.bam
stringtie -p 40 sort.bam --mix -t -f 0.1 -m 150 >stringtie.gtf

gtf_genome_to_cdna_fasta.pl stringtie.gtf genome.fa  >transcripts.fasta
gtf_to_alignment_gff3.pl stringtie.gtf > transcripts.gff3
TransDecoder.LongOrfs -t transcripts.fasta
makeblastdb -in db.fasta -dbtype prot -out dbname
blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db /PUBLIC/chicken/ref/anno/BY/00.data/protein/uniprotkb_organism_name_Gallus_2024_05_04.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 40 > blastp.outfmt6
hmmscan --cpu 40 --domtblout pfam.domtblout /PUBLIC/software/Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

#miniprot
#conda install bioconda::miniprot
#github:https://github.com/lh3/miniprot

miniprot -t16 -d genome.mpi genome.fa
miniprot -G 50k --gff -Iut 20 ./genome.mpi uniprotkb_organism_name_Gallus_2024_05_04.fasta >Gallus.pep.gff
miniprot_GFF_2_EVM_GFF3.py Gallus.pep.gff >Gallus.pep.gff3

#liftoff
#conda install bioconda::liftoff
#github:https://github.com/agshumate/Liftoff

liftoff -cds -polish -copies -flank 0.4 -sc 0.85 -g 7B.gff3 -o liftoff.gff3 -p 60 -chroms chr ./genome.fa Gal7B.sm.fa