#ab inito prediction
#augustus
#conda install bioconda::augustus
#github:https://github.com/Gaius-Augustus/Augustus

pblat -threads=60  -minIdentity=92 genome.fa all.flcdn.fa all.psl
pslCDnaFilter -maxAligns=1 all.psl cdna.f.psl
cat cdna.f.psl | sort -n -k 16,16 | sort -s -k 14,14 >cdna.f1.psl
blat2hints.pl --in=cdna.f1.psl --out=hints.E.gff
augustus --species=chicken --hintsfile=hints.E.gff --extrinsicCfgFile=extrinsic.M.RM.E.W.cfg  genome.maskN.fa > augustu.gff
augustus_GFF3_to_EVM_GFF3.pl augustu.gff >augustu.gff3