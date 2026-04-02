#!/usr/bin/env python3
"""
vg_surject_extract_tss.py
Extract the TSS (1-based) of each query(transcript) on the backbone from surjected BAM
usage:
  python3 vg_surject_extract_tss.py transcripts_surjected.sorted.bam > tx_to_backbone_TSS.tsv
Output:
  transcript_id \t backbone_chr \t tss_1based \t strand \t mapq \t aln_ref_len \t aln_cov_pct \t flag_notes
"""
import pysam,sys
bamfile = sys.argv[1]
sam = pysam.AlignmentFile(bamfile)
out = sys.stdout
seen = {}
for r in sam.fetch(until_eof=True):
    tid = r.query_name
    if r.is_unmapped:
        out.write(f"{tid}\tNO_MAPPING\t.\t.\t0\t0\t0\tUNMAPPED\n")
        continue
    # reference contig and mapping pos
    ref = sam.get_reference_name(r.reference_id)
    mapq = r.mapping_quality
    # compute reference aligned length (M, =, X, D, N contribute)
    aln_ref_len = 0
    for (c,l) in r.cigartuples or []:
        # cigartuple codes: M=0, I=1, D=2, N=3, S=4, H=5, P=6, = 7, X=8
        if c in (0,2,3,7,8):
            aln_ref_len += l
    # get reference start (0-based)
    rstart = r.reference_start
    # compute TSS depending on strand
    if r.is_reverse:
        # for reverse align, transcript 5' corresponds to alignment end on reference
        tend0 = rstart + aln_ref_len - 1   # 0-based last aligned base
        tss_1based = tend0 + 1
        strand = '-'
    else:
        tss_1based = rstart + 1
        strand = '+'
    # coverage pct: need transcript length
    qlen = r.query_length if r.query_length is not None else 0
    aln_cov_pct = (aln_ref_len/qlen*100) if qlen>0 else 0
    out.write(f"{tid}\t{ref}\t{tss_1based}\t{strand}\t{mapq}\t{aln_ref_len}\t{aln_cov_pct:.1f}\tOK\n")
sam.close()
