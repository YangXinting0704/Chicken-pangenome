#transcriptome-based predication
#RNA-seq
#conda install dnachun::star
#github:https://github.com/alexdobin/STAR https://github.com/gpertea/stringtie

sample=genome
STAR  --limitGenomeGenerateRAM=300000000000  --runMode genomeGenerate --runThreadN 20 --genomeDir ./ --genomeFastaFiles genome.fa
STAR --outSAMstrandField intronMotif  --runThreadN 20 --outFilterMismatchNoverLmax 0.02  --runMode alignReads --genomeDir /PUBLIC/chicken/ref/anno/D/  --alignIntronMin 20 --alignIntronMax 300000  --readFilesIn *_1.fq.gz *_2.fq.gz  --readFilesCommand zcat  --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outFileNamePrefix $sample
samtools merge -@ 70 merged.bam ./*/*bam
stringtie  -m 150 -t -f 0.1 -p 70 -o stringtie.gtf merged.bam

gtf_genome_to_cdna_fasta.pl stringtie.gtf genome.fa  >transcripts.fasta
gtf_to_alignment_gff3.pl stringtie.gtf > transcripts.gff3
TransDecoder.LongOrfs -t transcripts.fasta
makeblastdb -in db.fasta -dbtype prot -out dbname
blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db /PUBLIC/chicken/ref/anno/BY/00.data/protein/uniprotkb_organism_name_Gallus_2024_05_04.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 40 > blastp.outfmt6
hmmscan --cpu 40 --domtblout pfam.domtblout /PUBLIC/software/Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3