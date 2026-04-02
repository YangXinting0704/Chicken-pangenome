#!/usr/bin/env python3
# normalize_transcripts_to_pan.py
# 用法: python3 normalize_transcripts_to_pan.py transcripts_to_pan.gtf > transcripts_to_pan.normalized.gtf

import sys
import re

if len(sys.argv) < 2:
    print("Usage: python3 normalize_transcripts_to_pan.py transcripts_to_pan.gtf > normalized.gtf", file=sys.stderr)
    sys.exit(1)

fin = open(sys.argv[1])
out_lines = []
warn_lines = []
line_no = 0

# regex: try to capture final ::ChrName:start-end(strand)
# e.g. ....::Chr1:71883-74409(-)
pat = re.compile(r'(::)?(?P<pan>Chr[^\s:;]+):(?P<pstart>\d+)-(?P<pend>\d+)\((?P<pstrand>[-+])\)\s*$')
# We also accept lowercase chr or other names by relaxing: (?i:chr)
pat2 = re.compile(r'(::)?(?P<pan>[^:]+):(?P<pstart>\d+)-(?P<pend>\d+)\((?P<pstrand>[-+])\)\s*$', re.IGNORECASE)

for line in fin:
    line_no += 1
    if line.startswith('#') or line.strip() == '':
        sys.stdout.write(line)
        continue
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 9:
        # not a standard gtf line — pass through and warn
        warn_lines.append((line_no, "not 9 columns"))
        sys.stdout.write(line)
        continue
    seqname = fields[0]
    start = fields[3]
    end = fields[4]
    strand = fields[6]
    attrs = fields[8]

    # Try to extract last pan coord
    m = pat.search(seqname)
    if not m:
        m = pat2.search(seqname)
    if not m:
        # cannot parse, keep as is but log
        warn_lines.append((line_no, "cannot parse pan coord from seqname: " + seqname))
        sys.stdout.write(line)
        continue

    pan_chr = m.group('pan')
    pan_start = int(m.group('pstart'))
    pan_end = int(m.group('pend'))
    pan_strand = m.group('pstrand')

    # start/end in file may be relative offsets (small numbers). convert to absolute:
    try:
        rstart = int(start)
        rend = int(end)
    except:
        warn_lines.append((line_no, "start/end not ints"))
        sys.stdout.write(line)
        continue

    abs_start = pan_start + rstart - 1
    abs_end = pan_start + rend - 1

    # prefer pan_strand if present
    new_strand = pan_strand if pan_strand in ['+','-'] else strand

    # add attributes: orig_source and orig_exon_rel (to record original seqname and relative coords)
    # ensure attrs ends with ;
    if not attrs.strip().endswith(';'):
        attrs = attrs.strip() + ';'
    attrs = attrs + f' orig_source "{seqname}"; orig_exon_rel "{rstart}-{rend}";'

    # optionally: if original transcript_id lacks sample prefix and the seqname contains a prefix (like GalB41G000200::...), you could add orig_source to transcript id — we keep transcript_id as-is.
    fields[0] = pan_chr
    fields[3] = str(abs_start)
    fields[4] = str(abs_end)
    fields[6] = new_strand
    fields[8] = attrs

    sys.stdout.write('\t'.join(fields) + '\n')

# write warnings to stderr
for ln, msg in warn_lines:
    print(f"Warning: line {ln}: {msg}", file=sys.stderr)

fin.close()
