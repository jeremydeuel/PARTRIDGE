# script to extract sequences
# requires pandas and pysam
import pysam
import pandas as pd
import os
import gzip
import sys
from revcomp import revcomp
from qualityseq import QualitySeq
from consensus import find_consensus
TABLE_OF_EVIDENCE_READS = 'evidence.reads.txt.gz'
PATH_TO_BAMS = 'resources'
BAM_SUFFIX = '.bycoords.bam'
OUT_FILE = 'sequences.fq.gz'

if __name__ == '__main__':
    TABLE_OF_EVIDENCE_READS, \
    PATH_TO_BAMS, \
    OUT_FILE, \
    mouse = sys.argv[1:]
    mouse = int(mouse)
    d = pd.read_table(TABLE_OF_EVIDENCE_READS,sep=" ")

    print(f"processing mouse {mouse}")
    bam_path = os.path.join(PATH_TO_BAMS, f"{mouse}{BAM_SUFFIX}" )
    mouse_reads = d.loc[d['mouse']==mouse]
    reads = {}
    if not os.path.exists(bam_path):
        print(f"ERROR: bam file for mouse {mouse} not found at path {bam_path}")
        exit(0)
    rnames = set(mouse_reads.rname) #shortcut to speed things up
    with pysam.AlignmentFile(bam_path) as bam:
        for line in bam:
            if line.is_qcfail: continue
            if line.is_duplicate: continue

            if line.query_name in rnames:
                iap_alignments = mouse_reads.loc[mouse_reads.rname==line.query_name]
                iap_alignments = iap_alignments.loc[iap_alignments['pair.read'] == (1 if line.is_read1 else 2)]
                if not iap_alignments.empty:
                    seqname = iap_alignments.seqname.iloc[0]
                    end = iap_alignments['ltr.end'].iloc[0]
                    breakpoint = iap_alignments['breakpoint'].iloc[0]
                    read = iap_alignments['read'].iloc[0]
                    strand = iap_alignments['strand'].iloc[0]
                    clip = iap_alignments['clip'].iloc[0]
                    key = f"{seqname}:{breakpoint}:{end}"
                    if not key in reads:
                        reads[key] = {'breakpoints': [], 'mates': []}
                    if read == 'breakpoint':
                        if not line.reference_name == seqname:
                            continue #ignore this
                        if clip == "right":
                            offset, position = line.get_aligned_pairs(matches_only=True)[-1]
                            offset = offset + position - breakpoint
                            if offset < 0 or offset >= len(line.seq): continue
                            reads[key]['breakpoints'].append(QualitySeq(line.seq[offset:], line.query_qualities[offset:]))
                        else:
                            offset, position = line.get_aligned_pairs(matches_only=True)[0]
                            offset = offset + position - breakpoint
                            if offset < 0 or offset >= len(line.seq): continue
                            reads[key]['breakpoints'].append(QualitySeq(line.seq[:offset], line.query_qualities[:offset]).revcomp())

                    elif read == 'mate':
                        if line.is_read1 ^ line.is_reverse ^ (strand == '+') ^ (clip == 'right'):
                                reads[key]['mates'].append(QualitySeq(line.seq, line.query_qualities))
                        else:
                                reads[key]['mates'].append(QualitySeq(line.seq, line.query_qualities).revcomp())
    with (gzip.open(OUT_FILE, 'wt') as out_file):
        file_output = {}
        for breakpoint in reads.keys():
            breakpoint_consensus = find_consensus(reads[breakpoint]['breakpoints'])
            if len(breakpoint_consensus) < 8:
                print(f"WARNING: Found no consensus for breakpoint {breakpoint}")
                continue
            mates = reads[breakpoint]['mates']
            while len(mates):
                search_string = breakpoint_consensus[-8:]
                has_seq, new_mates = [], []
                for mate in mates:
                    if mate.has_once(search_string):
                        has_seq.append(mate)
                    else:
                        new_mates.append(mate)
                cons_seq = []
                mates = new_mates
                for seq in has_seq:
                    cons_seq.append(seq[(seq.find(search_string)+len(search_string)):])
                if len(cons_seq):
                    cons_seq = find_consensus(cons_seq)
                    if not len(cons_seq): break
                    #print(f"extended consensus: {breakpoint_consensus} | {cons_seq}")
                    breakpoint_consensus += cons_seq
                else:
                    break
            out_file.write(breakpoint_consensus.fastq(f"{breakpoint}:MOUSE{mouse}:CLIP"))
            i=0
            for mate in mates:
                i += 1
                out_file.write(mate.fastq(f"{breakpoint}:MOUSE{mouse}:MATE{i}"))




