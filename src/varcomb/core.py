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

    def remove_duplicates(self):
        return VCF(sorted(list(set(self.rows))))
