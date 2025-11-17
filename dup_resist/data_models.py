"""Data models used throughout Dup-Resist."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Optional


@dataclass
class GenePair:
    """Representation of a pair of collinear genes."""

    anchor1: str
    anchor2: str
    block_id: str
    pair_index: int
    score: Optional[float] = None
    e_value: Optional[float] = None
    orientation: Optional[str] = None
    ka: Optional[float] = None
    ks: Optional[float] = None

    def as_tuple(self) -> tuple[str, str]:
        return (self.anchor1, self.anchor2)


@dataclass
class CollinearityBlock:
    """A collinear block reported by MCScanX."""

    block_id: str
    score: Optional[float]
    e_value: Optional[float]
    chromosome1: Optional[str]
    chromosome2: Optional[str]
    pairs: List[GenePair] = field(default_factory=list)

    def add_pair(self, pair: GenePair) -> None:
        self.pairs.append(pair)

    def __iter__(self) -> Iterable[GenePair]:
        return iter(self.pairs)


@dataclass
class GeneConversionEvent:
    """Summary of a putative gene conversion event."""

    species: str
    block_id: str
    pair_index: int
    gene1: str
    gene2: str
    ks: float
    ka: Optional[float]
    block_median_ks: float
    block_zscore: Optional[float]
    detection_reason: str

    def to_row(self) -> List[str]:
        return [
            self.species,
            self.block_id,
            str(self.pair_index),
            self.gene1,
            self.gene2,
            f"{self.ks:.6f}",
            "" if self.ka is None else f"{self.ka:.6f}",
            f"{self.block_median_ks:.6f}",
            "" if self.block_zscore is None else f"{self.block_zscore:.3f}",
            self.detection_reason,
        ]
