#!/bin/bash 

set -e 

sample=$1
r1=$2
r2=$3

gtf="/ossfs/database/rna/index/V03_ChickenPan.gtf"
staridx="/ossfs/database/rna/index/"
refsalmonidx="/ossfs/database/rna/salmon/ChickenPan.salmon/"

#index
#V03_ChickenPan.fa
#V03_ChickenPan.gtf

threads=$(nproc)

cd /data
docker load < /ossfs/docker_images/star.tar.gz
#docker load < /ossfs/docker_images/salmon.tar.gz 

docker run -v /ossfs:/ossfs -v $(pwd):/data -w /data star:2.7.1a \
     STAR \
                --genomeDir ${staridx} \
                --sjdbGTFfile ${gtf} \
                --twopassMode Basic \
                --readFilesIn $r1 $r2 \
                --outFileNamePrefix /data/${sample}_ \
                --runThreadN ${threads} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterMismatchNoverLmax 0.04 \
                --outSAMunmapped Within \
                --chimSegmentMin 10 \
                --chimOutType Junctions \
                --chimOutJunctionFormat 1 \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 \
                --outSAMattributes NH HI AS nM NM ch \
                --alignSJDBoverhangMin 1 \
                --outSAMattrRGline ID:$sample SM:$sample

#samtools view -h -q 255 /data/${sample}_Aligned.sortedByCoord.out.bam | grep -v "vW:i:[2-7]" | samtools view -b > /data/${sample}_Aligned.sortedByCoord.out.filtered.bam

# 1. Expression quantification for exon (featureCounts v***)
docker run -v /ossfs:/ossfs -v /data/:/data -w /data localhost/PGCs \
    featureCounts -T ${threads} -p -t exon -g gene_id \
        -a ${gtf} \
        -o ${sample}_gene.tsv \
        /data/${sample}_Aligned.sortedByCoord.out.bam


# 2. Quantification for splicing events
cp /ossfs/soft/leafcutter-0.2.7.tar /opt/
cd /opt/ 
tar -xvf /opt/leafcutter-0.2.7.tar
export PATH=/opt/leafcutter-0.2.7/scripts/:$PATH
cd /data

# Step 1. BAM to junc
/opt/leafcutter-0.2.7/scripts/bam2junc.sh /data/${sample}_Aligned.sortedByCoord.out.bam  /data/${sample}_Aligned.sortedByCoord.out.bam.junc

# Step 2. Generating intron excision ratios with LeafCutter
python2 leafcutter_cluster.py -j juncfiles.txt -m 50 -p 0.001 -o tissue -l 500000


