###################################################################Extract transcript sequences
gffread ChickenPan.gff3 -g ChickenPan.fa -w ChickenPan_PanBaRT_transcripts.fa



###################################################################Projection of transcript sequences onto the pangenome graph
vg giraffe -Z pangenome.gbz -z pangenome.shortread.zipcodes -m pangenome.shortread.withzip.min -d pangenome.dist -f ChickenPan_PanBaRT_transcripts.fa > ChickenPan_PanBaRT_transcripts.gam
vg surject -x pangenome.pg -b ChickenPan_PanBaRT_transcripts.gam > ChickenPan_PanBaRT_transcripts.bam
python3 vg_surject_extract_tss.py ChickenPan_PanBaRT_transcripts.bam > ChickenPan_PanBaRT_transcripts.tsv
