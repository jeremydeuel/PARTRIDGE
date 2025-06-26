import pysam, os
import pandas as pd
from insertion import Insertion, Alignment, Multimap
BAM_DIR = "resources"
SOMATIC_INSERITONS = "somatic.insertions.txt.gz"
if __name__ == '__main__':
    sources = {}
    ins = pd.read_table(SOMATIC_INSERITONS, sep=" ")
    for f in os.listdir(BAM_DIR):
        if ".bam" in f:
            print(f)
            multimaps = {}
            with pysam.AlignmentFile(os.path.join(BAM_DIR, f)) as bam:
                for line in bam:
                    if line.is_unmapped: continue
                    seqname, pos, end, mouse, type = line.query_name.split(":")
                    id = ":".join((seqname, pos, end))
                    if id not in multimaps:
                        multimaps[id] = (Multimap(), {})

                    if type == "CLIP":
                        if not (end in ("3F","5R")) ^ line.is_forward:
                            alignment = Alignment(line.reference_name, line.reference_end, line.reference_start,
                                                not line.is_forward)
                        else:
                            alignment = Alignment(line.reference_name, line.reference_start, line.reference_end, line.is_forward)

                        multimaps[id][0].append(alignment)

                    else: #mate
                        if (end in ("3F", "5R")) ^ line.is_forward:
                            alignment = Alignment(line.reference_name, line.reference_end, line.reference_start,
                                            not line.is_forward)
                        else:
                            alignment = Alignment(line.reference_name, line.reference_start, line.reference_end,
                                                line.is_forward)
                        if not type in multimaps[id][1]:
                            multimaps[id][1][type] = Multimap([alignment])
                        else:
                            multimaps[id][1][type].append(alignment)

            for key in multimaps.keys():
                breakpoints = multimaps[key][0].collapse_breakpoints()
                if not len(breakpoints):
                    print(f"nothing left to match...")
                    break
                print(f"{key} with {len(breakpoints)} mappings for clipped part.")
                for type in multimaps[key][1].keys():
                    if len(multimaps[key][1][type]) < 10_000:
                        breakpoints = breakpoints.intersect(multimaps[key][1][type])
                #breakpoints = multimaps[key][0].collapse_breakpoints()
                multimaps[key] = breakpoints
    exit(0)
    for idx, row in ins.iterrows():
        if row.direction == "R":
            left_key = f"{row.seqname}:{row.LBP}:5R"
            right_key = f"{row.seqname}:{row.RBP}:3R"
        else:
            left_key = f"{row.seqname}:{row.LBP}:3F"
            right_key = f"{row.seqname}:{row.RBP}:5F"
        if left_key in sources.keys() and right_key in sources.keys():
            print(f"insertion {row.direction} {row.seqname}:{row.RBP}-{row.LBP}")
            left_matches = sources[left_key]
            right_matches = sources[right_key]
            matches = None
            for idx, l in left_matches['clip'].items():
                new_matches = set()
                for m in l:
                    if matches is None:
                        new_matches.append((m[0],m[1],m[2],m[3]))
                    else:
                        for n in matches:
                            if m[0]==n[0] and n[3]==m[3] and abs(m[1]-n[1])<10000:
                                new_matches.append((m[0],min(m[1],n[1]),max(m[2],n[2]),m[3]))
                if len(new_matches): matches = new_matches
            left_matches['clip'] = matches
            for idx, l in left_matches['clip'].items():
                new_matches = set()
                for m in l:
                    if matches is None:
                        new_matches.append((m[0],m[1],m[2],m[3]))
                    else:
                        for n in matches:
                            if m[0]==n[0] and n[3]==m[3] and abs(m[1]-n[1])<10000:
                                new_matches.append((m[0],min(m[1],n[1]),max(m[2],n[2]),m[3]))
                if len(new_matches): matches = new_matches

            for m in left_matches['clip']:
                for n in right_matches['clip']:
                    if m[0] == n[0] and m[3] != n[3]:
                        if abs(m[1]-n[1])<10000:
                            matches.append((n[0],min(m[1],n[1]), max(m[2],n[2]), '+' if n[3] else '-'))
            print(f"  .. have {len(matches)} after analysing breakpoints")
            # intersect with left mates
            for k, l in left_matches['mates'].items():
                new_matches = set()
                for m in l:
                    for n in matches:
                        if n[0] == m[0] and abs(n[1]-m[1]) < 10000:
                            new_matches.add((n[0],min(m[1],n[1]),max(m[2],n[2]), n[3]))
                if len(new_matches): matches = new_matches
            print(f"  .. have {len(matches)} after analysing left mates")
            # intersect with right mates
            for k, l in right_matches['mates'].items():
                new_matches = set()
                for m in l:
                    for n in matches:
                        if n[0] == m[0] and abs(n[1]-m[1]) < 10000:
                            new_matches.add((n[0],min(m[1],n[1]),max(m[2],n[2]), n[3]))
                if len(new_matches): matches = new_matches
            #intersect matches with eachother
            print(f"  .. have {len(matches)} after analysing right mates")
            new_matches = set()
            for i in range(len(matches)):
                m = matches[i]
                for j in range(i+1,len(matches)):
                    n = matches[j]
                    if n[0] == m[0] and abs(n[1]-m[1]) < 10000 and m[3]==n[3]:
                        nm = (n[0], min(m[1], n[1]), max(m[2], n[2]), n[3])
                        new_matches.add(nm)
            matches = list(new_matches)



            print(f"  found {len(matches)} matches: {matches}")
