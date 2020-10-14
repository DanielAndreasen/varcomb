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


class Info:
    def __init__(self, info):
        self.info = info

    def _convert_to_dict(self, info):
        if isinstance(info, str):
            self.info = dict()
            if info != '':
                for k in info.split(';'):
                    if '=' in k:
                        key, value = k.split('=')
                        self.info[key] = value
                    else:
                        self.info[k] = True

    def __eq__(self, other):
        return self.info == other.info

    def __len__(self):
        if isinstance(self.info, str):
            if ';' in self.info:
                return self.info.count(';') + 1
            return 0
        return len(self.info)

    def __getitem__(self, key):
        self._convert_to_dict(self.info)
        return self.info[key]

    def __setitem__(self, key, value):
        if isinstance(self.info, str):
            if ';' in self.info:
                self.info += f';{key}={value}'
            else:
                self.info += f'{key}={value}'
        else:
            self.info[key] = value

    def __str__(self):
        if isinstance(self.info, str):
            return self.info
        else:
            return ';'.join([f'{k}={v}' for k, v in self.info.items()])

    def keys(self):
        self._convert_to_dict(self.info)
        return self.info.keys()

    def values(self):
        self._convert_to_dict(self.info)
        return self.info.values()


@dataclass
class VCFrow:
    loc: Location
    id: str
    ref: str
    alt: str
    qual: str
    filter: str
    info: Info
    format: str
    samples: List[str]

    def __lt__(self, other) -> bool:
        return self.loc < other.loc

    def __hash__(self):
        return hash((self.loc, self.id, self.ref, self.alt, self.qual, self.filter, self.format))

    def _format_row(self):
        res = [self.loc.chrom, self.loc.pos, self.id, self.ref, self.alt, self.qual, self.filter, str(self.info), self.format]
        return '\t'.join(map(str, res + self.samples))


@dataclass
class VCF:
    rows: List[VCFrow]
    header: Optional[List[str]] = None

    def __str__(self):  # pragma: no cover
        o = ''
        for row in self.rows:
            o += row._format_row() + '\n'
        return o

    def __len__(self):
        return len(self.rows)

    def __add__(self, other):
        # TODO: Merge headers
        return VCF(self.rows + other.rows, header=self.header)

    def __getitem__(self, x: int) -> VCFrow:
        return self.rows[x]

    def __eq__(self, other):
        for row in self.rows:
            if row not in other.rows:
                return False
        return True

    def get_from_chrom(self, chrom):
        return VCF([row for row in self.rows if row.loc.chrom == chrom], header=self.header)

    def get_near_location(self, chrom, pos, tol=50):
        rows = []
        for row in self.get_from_chrom(chrom):
            if (pos - tol < row.loc.pos) and (row.loc.pos < pos + tol):
                rows.append(row)
        return VCF(rows, header=self.header)

    def remove_true_duplicates(self):
        return VCF(sorted(list(set(self.rows))), header=self.header)

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
        return VCF(sorted(list(set(rows))), header=self.header)

    def to_file(self, fname):
        header = '\n'.join(self.header) if self.header is not None else ''
        rows = []
        for row in self.rows:
            rows.append(row._format_row())

        rows = '\n'.join(rows)
        data = f'{header}\n{rows}'
        with open(fname.replace('.gz', ''), 'w') as f:
            f.write(data)

    def annotate(self, value: str):
        rows: List[VCFrow] = [row for row in self.rows]
        for i, row in enumerate(rows):
            row.info['Annotation'] = value
            rows[i] = row
        return VCF(sorted(rows), header=self.header)
