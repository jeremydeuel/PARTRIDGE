from typing import List

MAX_DISTANCE = 10_000

class Insertion:
    def __init__(self, seqname: str, LBP: int, RBP: int, is_forward: bool):
        self.seqname = seqname
        self.LBP = LBP
        self.RBP = RBP
        self.is_forward = is_forward


class Alignment:

    def __init__(self, seqname: str, breakpoint_start: int, alignment_end: int, is_forward: bool):
        self.seqname = seqname
        self.breakpoint_start = breakpoint_start
        self.alignment_end = alignment_end
        self.is_forward = is_forward

    def __str__(self) -> str:
        return f"alignment for {self.seqname}:{self.breakpoint_start}->{self.alignment_end} {'+' if self.is_forward else '-'}"

    def is_collapsible(self, other: 'Alignment') -> bool:
        if self.is_forward != other.is_forward:
            return False
        if self.seqname != other.seqname:
            return False
        if self.direction_rtl != other.direction_rtl:
            return False
        return abs(self.breakpoint_start - other.breakpoint_start) < 500

    def collapse(self, other: 'Alignment') -> 'Alignment':
        assert self.is_forward == other.is_forward
        assert self.seqname == other.seqname
        #assert self.direction_rtl == other.direction_rtl
        if self.direction_rtl: #reverse
            return Alignment(self.seqname, max(self.breakpoint_start, other.breakpoint_start), min(self.alignment_end, other.alignment_end), self.is_forward)
        else:
            return Alignment(self.seqname, min(self.breakpoint_start, other.breakpoint_start), max(self.alignment_end, other.alignment_end), self.is_forward)

    @property
    def direction_rtl(self):
        """alignment is left to right direction"""
        return self.breakpoint_start > self.alignment_end



class Multimap:

    def __init__(self, alignment: List[Alignment] = None):
        if alignment is None:
            self.alignment = []
        else:
            self.alignment = alignment

    def __len__(self) -> int:
        return len(self.alignment)

    def append(self, alignment: Alignment):
        self.alignment.append(alignment)

    def collapse_breakpoints(self):
        if len(self.alignment) < 2:
            print(f"nothing to collapse, len={len(self.alignment)}")
            return self
        collapsed = Multimap()
        excluded = set()
        for i in range(len(self.alignment)):
            x = self.alignment[i]
            for j in range(i+1, len(self.alignment)):
                if j in excluded: continue
                if x.is_collapsible(self.alignment[j]):
                    x = x.collapse(self.alignment[j])
                    excluded.add(j)
            collapsed.append(x)
        #print(f" collapsed {len(self.alignment)} maps to {len(collapsed)} maps.")
        return collapsed

    def intersect(self, other: 'Multimap') -> 'Multimap':
        intersection = Multimap()
        print(f"intersecting {len(self)} with {len(other)} multimaps, have {len(intersection)}")
        for i in range(len(self)):
            for j in range(len(other)):
                if self.alignment[i].seqname == other.alignment[j].seqname and \
                    self.alignment[i].is_forward == other.alignment[j].is_forward:
                    #print(f"{self.alignment[i]} <-> {other.alignment[j]}")
                    if self.alignment[i].direction_rtl:
                        if other.alignment[j].breakpoint_start <= self.alignment[i].breakpoint_start and \
                                self.alignment[i].breakpoint_start - other.alignment[j].breakpoint_start < 5_000:
                            intersection.append(Alignment(self.alignment[i].seqname, self.alignment[i].breakpoint_start, min(self.alignment[i].alignment_end, other.alignment[j].alignment_end, other.alignment[j].breakpoint_start), self.alignment[i].is_forward))
                    else:
                        if other.alignment[j].breakpoint_start >= self.alignment[i].breakpoint_start and \
                            other.alignment[j].breakpoint_start - self.alignment[i].breakpoint_start < 5_000:
                            intersection.append(Alignment(self.alignment[i].seqname, self.alignment[i].breakpoint_start, max(self.alignment[i].alignment_end, other.alignment[j].alignment_end, other.alignment[j].breakpoint_start), self.alignment[i].is_forward))

        #
        print(f".. reduced to {len(intersection)}")
        return intersection

