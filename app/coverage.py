"""Paired-end coverage check against a small reference sequence.

For a succeeded pairing task this streams the published R1/R2 result files
(read-only), places every read pair on a single ACGT reference (<= 10,000
bp) by exact matching and accumulates per-base depth.

A pair contributes depth only when it has *exactly one* compliant
placement: the left end (R1) matches the reference in forward orientation,
the right end (R2) matches as reverse complement, start positions are
strictly increasing and the outermost covered span is at most 500 bases.
Each end then adds +1 to every base it covers, so an overlap of the two
ends is counted twice. Pairs with no compliant placement and pairs with
several compliant placements are counted separately and add no depth.

The report is the list of contiguous below-threshold intervals on the
reference as 1-based closed intervals plus the three pair counters.

Performance: match positions are found with ``bytes.find`` and cached per
distinct read sequence inside a fixed budget of int32 positions, and the
placement search stops after two hits, so repetitive references with
astronomically many candidate combinations are classified in milliseconds
without ever materializing the combinations. The working set is a handful
of small arrays: the reference (<= 10 KB), one difference array
(<= 10,001 ints) and at most ``COVERAGE_CACHE_POSITIONS`` cached int32
positions (default 1M -> ~4 MiB), all far below the 256 MiB cap.
"""

from __future__ import annotations

import bisect
import itertools
import re
from array import array
from pathlib import Path
from typing import Iterator

from . import config
from .fastq import FastqFormatError, _iter_raw_records, extract_read_id


class CoverageInputError(ValueError):
    """Raised when any coverage-check input violates the contract."""


_RC_TABLE = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
_SEQUENCE_LINE = re.compile(r"^[ACGTacgt]+$")


def reverse_complement(seq: bytes) -> bytes:
    """Reverse complement of an ACGT sequence (other bytes pass through)."""
    return seq.translate(_RC_TABLE)[::-1]


def parse_reference_fasta(text: str) -> tuple[str, str]:
    """Validate a single-sequence ACGT FASTA; return (name, sequence).

    Lowercase bases are accepted and normalized to upper case. Any
    deviation (missing header, a second sequence, non-ACGT bases, empty or
    over-long sequence) raises CoverageInputError.
    """
    if not text or not text.strip():
        raise CoverageInputError("reference FASTA is empty")
    lines = text.splitlines()
    first = next((i for i, line in enumerate(lines) if line.strip()), None)
    if first is None or not lines[first].lstrip().startswith(">"):
        raise CoverageInputError("reference FASTA must start with a '>' header line")
    header = lines[first].lstrip()[1:].strip()
    name = header.split(None, 1)[0] if header else ""

    chunks: list[str] = []
    total = 0
    for lineno, raw in enumerate(lines[first + 1:], start=first + 2):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            raise CoverageInputError(
                "reference FASTA must hold exactly one sequence; "
                f"second header found at line {lineno}"
            )
        if not _SEQUENCE_LINE.match(line):
            raise CoverageInputError(
                f"reference sequence line {lineno} contains non-ACGT characters"
            )
        chunks.append(line.upper())
        total += len(line)
        if total > config.COVERAGE_MAX_REFERENCE_BASES:
            raise CoverageInputError(
                f"reference sequence exceeds the "
                f"{config.COVERAGE_MAX_REFERENCE_BASES} bp limit"
            )
    if total == 0:
        raise CoverageInputError("reference FASTA contains no sequence")
    return name, "".join(chunks)


class _MatchIndex:
    """All start positions of a read within the reference, cached per read.

    Positions are stored as sorted int32 arrays. The cache is bounded by a
    total number of cached positions, so even thousands of distinct highly
    repetitive reads cannot exhaust memory; once the budget is spent the
    index simply recomputes without caching.
    """

    __slots__ = ("_ref", "_cache", "_cached", "_budget")

    def __init__(self, ref: bytes) -> None:
        self._ref = ref
        self._cache: dict[bytes, array] = {}
        self._cached = 0
        self._budget = config.COVERAGE_CACHE_POSITIONS

    def positions(self, seq: bytes) -> array:
        cached = self._cache.get(seq)
        if cached is not None:
            return cached
        found = array("i")
        find = self._ref.find
        start = 0
        while True:
            pos = find(seq, start)
            if pos < 0:
                break
            found.append(pos)
            start = pos + 1
        new_total = self._cached + len(found)
        if new_total <= self._budget:
            self._cache[seq] = found
            self._cached = new_total
        return found


def _placements_capped(
    p1_positions: array,
    p2_positions: array,
    len2: int,
    cap: int = 2,
) -> list[tuple[int, int]]:
    """Return up to ``cap`` compliant (p1, p2) placements of one read pair.

    A compliant placement has p1 < p2 (strictly increasing starts) and an
    outermost covered span p2 + len2 - p1 <= COVERAGE_MAX_SPAN. Callers only
    ever need to distinguish 0 / 1 / many, so the search stops at ``cap``.
    """
    window = config.COVERAGE_MAX_SPAN - len2
    out: list[tuple[int, int]] = []
    if window < 1 or not p1_positions or not p2_positions:
        return out
    for p1 in p1_positions:
        lo = bisect.bisect_right(p2_positions, p1)
        hi = bisect.bisect_right(p2_positions, p1 + window)
        for k in range(lo, hi):
            out.append((p1, p2_positions[k]))
            if len(out) >= cap:
                return out
    return out


def _iter_pairs(r1_path: Path, r2_path: Path) -> Iterator[tuple[bytes, bytes]]:
    """Yield (seq1, seq2) per pair while validating both FASTQ streams.

    Structure, read-id pairing and per-end length are checked exactly; any
    violation raises CoverageInputError so the whole check fails instead of
    returning partial intervals.
    """
    sentinel = object()
    streams = itertools.zip_longest(
        _iter_raw_records(r1_path, "R1"),
        _iter_raw_records(r2_path, "R2"),
        fillvalue=sentinel,
    )
    pair_no = 0
    for rec1, rec2 in streams:
        if rec1 is sentinel or rec2 is sentinel:
            raise CoverageInputError(
                "R1 and R2 result files hold different numbers of records"
            )
        pair_no += 1
        if pair_no > config.COVERAGE_MAX_PAIRS:
            raise CoverageInputError(
                f"more than {config.COVERAGE_MAX_PAIRS} read pairs; "
                "the coverage check processes at most "
                f"{config.COVERAGE_MAX_PAIRS} pairs"
            )
        id1 = extract_read_id(rec1[0].rstrip(b"\r\n"))
        id2 = extract_read_id(rec2[0].rstrip(b"\r\n"))
        if id1 != id2:
            raise CoverageInputError(
                f"pair {pair_no}: read id mismatch '{id1}' vs '{id2}'"
            )
        seq1 = rec1[1].rstrip(b"\r\n")
        seq2 = rec2[1].rstrip(b"\r\n")
        for label, seq in (("R1", seq1), ("R2", seq2)):
            if not (
                config.COVERAGE_MIN_READ_BASES
                <= len(seq)
                <= config.COVERAGE_MAX_READ_BASES
            ):
                raise CoverageInputError(
                    f"pair {pair_no}: {label} read length {len(seq)} bp is "
                    f"outside the {config.COVERAGE_MIN_READ_BASES}-"
                    f"{config.COVERAGE_MAX_READ_BASES} bp range"
                )
        yield seq1, seq2


def check_pair_coverage(
    r1_path: Path,
    r2_path: Path,
    reference_fasta: str,
    depth_threshold: int,
) -> dict:
    """Run the paired coverage check and return the report as a dict.

    Raises CoverageInputError on any invalid input (reference, threshold,
    FASTQ structure, read ids, lengths, pair count). The result files are
    only read, never modified.
    """
    name, reference = parse_reference_fasta(reference_fasta)
    if not (
        config.COVERAGE_MIN_DEPTH_THRESHOLD
        <= depth_threshold
        <= config.COVERAGE_MAX_DEPTH_THRESHOLD
    ):
        raise CoverageInputError(
            f"depth threshold must be within "
            f"{config.COVERAGE_MIN_DEPTH_THRESHOLD}-"
            f"{config.COVERAGE_MAX_DEPTH_THRESHOLD}"
        )

    ref = reference.encode("ascii")
    ref_len = len(ref)
    index = _MatchIndex(ref)
    diff = [0] * (ref_len + 1)
    placed = no_placement = multiple = 0

    try:
        for seq1, seq2 in _iter_pairs(r1_path, r2_path):
            p1_positions = index.positions(seq1)
            p2_positions = index.positions(reverse_complement(seq2))
            hits = _placements_capped(p1_positions, p2_positions, len(seq2))
            if not hits:
                no_placement += 1
            elif len(hits) > 1:
                multiple += 1
            else:
                placed += 1
                p1, p2 = hits[0]
                diff[p1] += 1
                diff[p1 + len(seq1)] -= 1
                diff[p2] += 1
                diff[p2 + len(seq2)] -= 1
    except FastqFormatError as exc:
        raise CoverageInputError(str(exc)) from exc

    intervals: list[dict] = []
    depth = 0
    run_start: int | None = None
    for i in range(ref_len):
        depth += diff[i]
        if depth < depth_threshold:
            if run_start is None:
                run_start = i
        elif run_start is not None:
            # 0-based run [run_start, i-1] -> 1-based closed interval.
            intervals.append({"start": run_start + 1, "end": i})
            run_start = None
    if run_start is not None:
        intervals.append({"start": run_start + 1, "end": ref_len})

    return {
        "reference_name": name,
        "reference_length": ref_len,
        "depth_threshold": depth_threshold,
        "pairs_total": placed + no_placement + multiple,
        "pairs_placed": placed,
        "pairs_no_placement": no_placement,
        "pairs_multiple_placements": multiple,
        "low_coverage_intervals": intervals,
    }
