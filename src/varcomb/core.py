from dataclasses import dataclass
from typing import List, Optional

from varcomb.exceptions import LocationShiftError


def as_int(x):
    try:
        xp = int(x)
    except ValueError:
        xp = x.upper()
    return xp


@dataclass
class Location:
    chrom: str
    pos: int

    def _trim_chrom(self):
        if self.chrom.lower().startswith('chr'):
            return as_int(self.chrom[3:])
        return as_int(self.chrom)

    def __sub__(self, other) -> Optional[int]:
        if self.chrom != other.chrom:
            return None
        return self.pos - other.pos

    def __add__(self, shift: int):
        if not isinstance(shift, int):
            raise LocationShiftError(f'Shifting is only allowed with an integer. "{shift}" was parsed')
        return Location(self.chrom, self.pos + shift)

    def __lt__(self, other):
        if self.chrom == other.chrom:
            return self.pos < other.pos
        chrom1 = self._trim_chrom()
        chrom2 = other._trim_chrom()
        if isinstance(chrom1, str) and isinstance(chrom2, int):
            return False
        elif isinstance(chrom1, int) and isinstance(chrom2, str):
            return True
        # Both are int or str
        return chrom1 < chrom2

    def __hash__(self):
        return hash(self.chrom) + hash(self.pos)

    def shift(self, shift: int):
        return self + shift


@dataclass
class VCFrow:
    loc: Location
    id: str
    ref: str
    alt: str
    qual: str
    filter: str
    info: str
    format: str
    samples: List[str]

    def __lt__(self, other) -> bool:
        return self.loc < other.loc

    def __hash__(self):
        return hash((self.loc, self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format))


@dataclass
class VCF:
    rows: List[VCFrow]

    def __len__(self):
        return len(self.rows)

    def __add__(self, other):
        return VCF(self.rows + other.rows)

    def __getitem__(self, x: int) -> VCFrow:
        return self.rows[x]

    def __eq__(self, other):
        for row in self.rows:
            if row not in other.rows:
                return False
        return True

    def get_from_chrom(self, chrom):
        return VCF([row for row in self.rows if row.loc.chrom == chrom])

    def get_near_location(self, chrom, pos, tol=50):
        rows = []
        for row in self.get_from_chrom(chrom):
            if (pos - tol < row.loc.pos) and (row.loc.pos < pos + tol):
                rows.append(row)
        return VCF(rows)

    def remove_true_duplicates(self):
        return VCF(sorted(list(set(self.rows))))

    def remove_loc_dup(self):
        rows: List[VCFrow] = [row for row in self.rows]
        completed = False
        while not completed:
            for row in rows:
                mask = [row.loc == r.loc for r in rows]
                if sum(mask) > 1:
                    mask[mask.index(True)] = False
                    rows = [row for row, m in zip(rows, mask) if not m]
                    break
                if sum(mask) == 1:
                    completed = True
            else:
                break
        return VCF(sorted(list(set(rows))))
