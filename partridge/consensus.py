# Geminate - Genomic events of mutation through insertional alterations by transposable elements
#
# Copyright (C) 2025 Jeremy Deuel <jeremy.deuel@usz.ch>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


from typing import Tuple, Iterable
from qualityseq import QualitySeq


def find_consensus(seqs: Iterable[QualitySeq]) -> Tuple[QualitySeq]:
    """
    Find consensus sequence of an iterable of sequences, that must be left-aligned
    returns a tuple of
        - list of integers: number of sequences matching consensus base at each position
        - list of floats: Fraction of sequences matching consensus sequences at each position
        - string: Consensus sequence
    """
    if len(seqs) == 0:
        return QualitySeq('', [])
    if len(seqs) == 1:
        return seqs[0]
    lengths = [len(s) for s in seqs]
    seqs = [s.upper() for s in seqs] #ensure upper case
    consensus_seq = ""
    consensus_score = []
    consensus_n = []
    for position in range(max(lengths)):
        base_stat = {
            'A': 0, 'T': 0, 'G': 0, 'C': 0
        }
        for sequence in range(len(lengths)):
            if lengths[sequence] <= position: continue
            base = seqs[sequence].seq(position)
            if not base in base_stat.keys():
                continue
            base_stat[base] += seqs[sequence].qual(position)
        n_sequences = sum(base_stat.values())
        n_sorted = sorted(base_stat.items(), key=lambda x: x[1], reverse=True)
        delta_first_two_hits = n_sorted[0][1]-n_sorted[1][1]
        if delta_first_two_hits > 25:
            consensus_score.append(min(40+int(delta_first_two_hits/37),delta_first_two_hits))
            consensus_seq += n_sorted[0][0]
        else:
            #consensus_score.append(0)
            #consensus_seq += "N"
            break #only extract seq to the first ambigous base.
            #consensus_score.append(0)
            #consensus_n.append(0)
            #consensus_seq += "N"
    return QualitySeq(consensus_seq, consensus_score)

