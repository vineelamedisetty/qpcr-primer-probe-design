#!/usr/bin/env python3
"""
qPCR primer + hydrolysis (TaqMan) probe design for pathogen detection.

Design flow
-----------
  FASTA  ->  per-record sequence
         ->  enumerate primer windows
         ->  filter (length, GC, nearest-neighbour Tm, runs, ambiguity,
                     3' GC clamp, hairpin, self-dimer)
         ->  pair forward x reverse within the amplicon window
             (start-indexed, so the pairing is linear in sequence length)
         ->  check cross-dimer, place a probe inside the amplicon interior
         ->  score, de-duplicate overlapping assays, rank, write CSV

Thermodynamics
--------------
Tm uses the SantaLucia (1998) unified nearest-neighbour parameters with the
Owczarzy (2004) monovalent-salt correction, at user-set primer and Na+
concentrations. The GC-content formula used by many teaching scripts
(64.9 + 41*(GC-16.4)/N) is sequence-order blind and ignores salt and strand
concentration; on a real target it misplaces Tm by several degrees.

Scope / limitation
------------------
This tool checks THERMODYNAMICS and STRUCTURE only. It does NOT check
SPECIFICITY. Before any wet-lab use, every candidate must be screened against
the host genome, near-neighbour organisms and the target's own strain
diversity (e.g. BLASTn vs nt, or Primer-BLAST). `--report-fasta` writes the
candidates out in FASTA form for exactly that purpose.

Dependencies: Python >= 3.8. pandas optional (falls back to csv module).
Usage: python qpcr_design.py                      (uses pathogen.fasta beside this file)
       python qpcr_design.py --fasta other.fasta --out assays.csv
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import sys
from dataclasses import dataclass, asdict
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# Folder this script lives in. Defaults are resolved against it, so the
# script finds pathogen.fasta no matter which directory it is run from
# (VS Code's Run button starts in the workspace root, not here).
HERE = Path(__file__).resolve().parent

# ============================================================================
# 1. DESIGN PARAMETERS  (edit here, or override on the command line)
# ============================================================================

# --- primers ---
MIN_PRIMER_LEN, MAX_PRIMER_LEN = 18, 25
MIN_PRIMER_GC, MAX_PRIMER_GC = 40.0, 60.0
MIN_PRIMER_TM, MAX_PRIMER_TM = 58.0, 62.0
MAX_PRIMER_TM_DIFF = 2.0          # forward vs reverse Tm mismatch

# --- probe (5' nuclease / TaqMan) ---
MIN_PROBE_LEN, MAX_PROBE_LEN = 18, 30
MIN_PROBE_GC, MAX_PROBE_GC = 30.0, 80.0
PROBE_TM_ABOVE_PRIMER = (6.0, 12.0)   # probe Tm must exceed the higher primer
                                      # Tm by this much (classic rule: ~8-10 C)
PROBE_PRIMER_GAP = 1              # min bases between a primer 3' end and probe

# --- amplicon ---
MIN_AMPLICON, MAX_AMPLICON = 70, 150
OPTIMAL_AMPLICON = 100

# --- structure ---
MAX_HOMOPOLYMER = 4               # reject 5+ identical bases in a row
MAX_SELF_DIMER = 8                # max complementary bases, whole molecule
MAX_3P_DIMER = 3                  # max complementary bases at the 3' end
MAX_HAIRPIN_STEM = 4              # max stem length with a >=3 nt loop
MIN_GC_CLAMP, MAX_GC_CLAMP = 1, 3  # G/C count in the last 5 bases

# --- reaction conditions used for Tm ---
NA_CONC = 0.05                    # monovalent cation, mol/L
PRIMER_CONC = 250e-9              # mol/L
PROBE_CONC = 250e-9               # mol/L

# --- reporting ---
TOP_RESULTS = 10
MIN_ASSAY_SEPARATION = 20         # bp between reported forward primer starts,
                                  # so the table is not one assay ten times

VALID_BASES = frozenset("ACGT")

# Full IUPAC complement: R<->Y, K<->M, B<->V, D<->H, S/W/N self-complementary.
IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYKMBVDHNSWacgtrykmbvdhnsw",
    "TGCAYRMKVBHDNSWtgcayrmkvbhdnsw",
)


# ============================================================================
# 2. THERMODYNAMICS
# ============================================================================

# SantaLucia (1998) unified parameters: dH kcal/mol, dS cal/(mol*K)
NN_PARAMS: Dict[str, Tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)    # helix initiation with a terminal G/C
INIT_AT = (2.3, 4.1)     # helix initiation with a terminal A/T
R_GAS = 1.9872           # cal/(mol*K)


def gc_percent(seq: str) -> float:
    """G+C as a percentage. Only unambiguous G/C are counted."""
    if not seq:
        return 0.0
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def tm_nearest_neighbour(seq: str,
                         strand_conc: float = PRIMER_CONC,
                         na_conc: float = NA_CONC) -> float:
    """Melting temperature in Celsius (SantaLucia 1998 + Owczarzy 2004 salt).

    Assumes a non-self-complementary oligo in excess over its template, which
    is the normal PCR case, hence the CT/4 term. Results are cached: the same
    window is evaluated many times while primer pairs are enumerated.
    """
    return _tm_cached(seq.upper(), strand_conc, na_conc)


@lru_cache(maxsize=262144)
def _tm_cached(seq: str, strand_conc: float, na_conc: float) -> float:
    if len(seq) < 2 or not set(seq) <= VALID_BASES:
        raise ValueError(f"Tm needs an unambiguous A/C/G/T oligo, got {seq!r}")

    d_h = d_s = 0.0
    for i in range(len(seq) - 1):
        h, s = NN_PARAMS[seq[i:i + 2]]
        d_h += h
        d_s += s
    for terminal in (seq[0], seq[-1]):
        h, s = INIT_GC if terminal in "GC" else INIT_AT
        d_h += h
        d_s += s

    tm_kelvin = (d_h * 1000.0) / (d_s + R_GAS * math.log(strand_conc / 4.0))

    # Owczarzy et al. 2004 monovalent salt correction (1/Tm form)
    f_gc = (seq.count("G") + seq.count("C")) / len(seq)
    ln_na = math.log(na_conc)
    inv_tm = (1.0 / tm_kelvin
              + (4.29 * f_gc - 3.95) * 1e-5 * ln_na
              + 9.40e-6 * ln_na ** 2)
    return 1.0 / inv_tm - 273.15


@lru_cache(maxsize=262144)
def reverse_complement(seq: str) -> str:
    """Reverse complement, IUPAC-aware (R->Y, Y->R, K->M, ...)."""
    return seq.upper().translate(IUPAC_COMPLEMENT)[::-1]


# ============================================================================
# 3. STRUCTURE CHECKS
# ============================================================================

_COMPLEMENT_PAIRS = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}


def _is_pair(a: str, b: str) -> bool:
    return (a, b) in _COMPLEMENT_PAIRS


@lru_cache(maxsize=262144)
def max_duplex_score(seq_a: str, seq_b: str) -> int:
    """Largest number of complementary bases over all ungapped alignments.

    seq_b is read 3'->5' against seq_a read 5'->3', which is how two oligos
    anneal. Used for self-dimer (a vs a) and cross-dimer (fwd vs rev).
    """
    b_rev = seq_b[::-1]
    best = 0
    for shift in range(-(len(seq_a) - 1), len(b_rev)):
        score = 0
        for i, base in enumerate(seq_a):
            j = i - shift
            if 0 <= j < len(b_rev) and _is_pair(base, b_rev[j]):
                score += 1
        best = max(best, score)
    return best


@lru_cache(maxsize=262144)
def three_prime_duplex_score(seq_a: str, seq_b: str, window: int = 5) -> int:
    """Complementarity involving the 3' terminus of seq_a.

    3' dimers are the dangerous ones: they are extendable by the polymerase
    and produce primer-dimer product that consumes reagent and, in SYBR
    formats, signal.
    """
    tail = seq_a[-window:]
    b_rev = seq_b[::-1]
    best = 0
    for shift in range(-(len(tail) - 1), len(b_rev)):
        run = 0
        for i, base in enumerate(tail):
            j = i - shift
            if 0 <= j < len(b_rev) and _is_pair(base, b_rev[j]):
                run += 1
            else:
                run = 0
            best = max(best, run)
    return best


@lru_cache(maxsize=262144)
def max_hairpin_stem(seq: str, min_loop: int = 3) -> int:
    """Longest self-complementary stem that can fold back with a >=3 nt loop."""
    best = 0
    n = len(seq)
    for i in range(n):
        for j in range(i + min_loop + 1, n):
            stem = 0
            while (i - stem >= 0 and j + stem < n
                   and (j + stem) - (i - stem) - 1 >= min_loop
                   and _is_pair(seq[i - stem], seq[j + stem])):
                stem += 1
            best = max(best, stem)
    return best


def has_long_homopolymer(seq: str, limit: int = MAX_HOMOPOLYMER) -> bool:
    """True if any base repeats more than `limit` times consecutively."""
    run, previous = 1, ""
    for base in seq:
        run = run + 1 if base == previous else 1
        if run > limit:
            return True
        previous = base
    return False


def gc_clamp_count(seq: str, window: int = 5) -> int:
    """G/C bases in the 3'-terminal window. 1-2 is the usual target."""
    tail = seq[-window:]
    return tail.count("G") + tail.count("C")


# ============================================================================
# 4. FASTA I/O
# ============================================================================

@dataclass
class FastaRecord:
    identifier: str
    description: str
    sequence: str


def read_fasta(path: str) -> List[FastaRecord]:
    """Parse a FASTA file into records.

    Every record is kept separate. Concatenating records into one string (a
    common shortcut) invents chimeric primers that span a junction present in
    no real molecule.
    """
    records: List[FastaRecord] = []
    header: Optional[str] = None
    chunks: List[str] = []

    def flush() -> None:
        if header is not None:
            ident, _, desc = header.partition(" ")
            records.append(FastaRecord(ident or "unnamed", desc.strip(),
                                       "".join(chunks).upper()))

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header, chunks = line[1:].strip(), []
            elif header is None:
                raise ValueError(f"{path}: sequence data before any '>' header")
            else:
                chunks.append(line)
    flush()

    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    for record in records:
        if not record.sequence:
            raise ValueError(f"{record.identifier}: empty sequence")
    return records


def summarise_record(record: FastaRecord) -> str:
    ambiguous = sum(1 for b in record.sequence if b not in VALID_BASES)
    note = f", {ambiguous} ambiguous base(s) excluded from design" if ambiguous else ""
    return f"{record.identifier}: {len(record.sequence)} bp{note}"


# ============================================================================
# 5. CANDIDATE GENERATION AND FILTERING
# ============================================================================

@dataclass
class Oligo:
    sequence: str      # as it appears on the plus strand
    start: int         # 0-based plus-strand start
    length: int
    gc: float
    tm: float

    @property
    def end(self) -> int:
        """Exclusive plus-strand end coordinate."""
        return self.start + self.length


def passes_primer_filters(seq: str, tm: float, gc: float) -> bool:
    """Composition, thermodynamic and structural gate for a single primer."""
    if not MIN_PRIMER_GC <= gc <= MAX_PRIMER_GC:
        return False
    if not MIN_PRIMER_TM <= tm <= MAX_PRIMER_TM:
        return False
    if has_long_homopolymer(seq):
        return False
    if not MIN_GC_CLAMP <= gc_clamp_count(seq) <= MAX_GC_CLAMP:
        return False
    if max_duplex_score(seq, seq) > MAX_SELF_DIMER:
        return False
    if three_prime_duplex_score(seq, seq) > MAX_3P_DIMER:
        return False
    if max_hairpin_stem(seq) > MAX_HAIRPIN_STEM:
        return False
    return True


def find_candidate_primers(sequence: str) -> List[Oligo]:
    """All plus-strand windows that survive the primer filters.

    Windows containing anything other than A/C/G/T are skipped: an N or an
    IUPAC code means the base is not known, so its Tm and GC cannot be
    computed and a primer placed there may not match the real target.
    """
    candidates: List[Oligo] = []
    n = len(sequence)
    for start in range(n - MIN_PRIMER_LEN + 1):
        for length in range(MIN_PRIMER_LEN, MAX_PRIMER_LEN + 1):
            end = start + length
            if end > n:
                break
            window = sequence[start:end]
            if not set(window) <= VALID_BASES:
                break  # an ambiguous base here poisons every longer window too
            gc = gc_percent(window)
            if not MIN_PRIMER_GC <= gc <= MAX_PRIMER_GC:
                continue
            tm = tm_nearest_neighbour(window, PRIMER_CONC, NA_CONC)
            if not MIN_PRIMER_TM <= tm <= MAX_PRIMER_TM:
                continue
            if passes_primer_filters(window, tm, gc):
                candidates.append(Oligo(window, start, length, gc, tm))
    return candidates


def index_by_start(primers: Sequence[Oligo], length: int) -> List[List[Oligo]]:
    """Bucket primers by plus-strand start position.

    This is what keeps pairing tractable. Comparing every primer with every
    other is O(m^2) and pointless, because only primers within MAX_AMPLICON
    of each other can ever form a product.
    """
    buckets: List[List[Oligo]] = [[] for _ in range(length + 1)]
    for primer in primers:
        buckets[primer.start].append(primer)
    return buckets


# ============================================================================
# 6. PROBE DESIGN
# ============================================================================

@dataclass
class Probe:
    sequence: str      # as synthesised (may be the minus-strand form)
    start: int         # 0-based plus-strand start of its binding site
    end: int           # exclusive
    tm: float
    gc: float
    c_minus_g: int     # C count minus G count; hydrolysis probes prefer > 0


def probe_structure_ok(seq: str, gc: float) -> bool:
    """Everything about a probe that does not depend on the primer pair."""
    if seq.startswith("G"):
        return False                      # 5' G quenches the reporter dye
    if not MIN_PROBE_GC <= gc <= MAX_PROBE_GC:
        return False
    if has_long_homopolymer(seq):
        return False
    if max_hairpin_stem(seq) > MAX_HAIRPIN_STEM:
        return False
    if max_duplex_score(seq, seq) > MAX_SELF_DIMER:
        return False
    return True


def build_probe_candidates(sequence: str) -> Tuple[List[Probe], List[int]]:
    """Every acceptable probe in the record, computed once.

    Only the Tm-relative-to-primers rule depends on the pair, so it is applied
    at selection time. Building this list once instead of rescanning inside
    every primer pair is what keeps the run linear in sequence length.
    Returns the candidates sorted by start plus their start positions, for
    binary search.
    """
    candidates: List[Probe] = []
    n = len(sequence)
    for start in range(n - MIN_PROBE_LEN + 1):
        for length in range(MIN_PROBE_LEN, MAX_PROBE_LEN + 1):
            end = start + length
            if end > n:
                break
            plus = sequence[start:end]
            if not set(plus) <= VALID_BASES:
                break  # an ambiguous base poisons every longer window too
            # The probe may sit on either strand.
            for oligo in (plus, reverse_complement(plus)):
                gc = gc_percent(oligo)
                if not MIN_PROBE_GC <= gc <= MAX_PROBE_GC:
                    continue
                if not probe_structure_ok(oligo, gc):
                    continue
                candidates.append(Probe(
                    sequence=oligo,
                    start=start,
                    end=end,
                    tm=tm_nearest_neighbour(oligo, PROBE_CONC, NA_CONC),
                    gc=gc,
                    c_minus_g=oligo.count("C") - oligo.count("G"),
                ))
    candidates.sort(key=lambda p: p.start)
    return candidates, [p.start for p in candidates]


def select_probe(candidates: Sequence[Probe], starts: Sequence[int],
                 forward: "Oligo", reverse: "Oligo") -> Optional[Probe]:
    """Best probe lying strictly between the two primer footprints.

    The probe must melt above both primers (classic rule ~8-10 C) so that it
    is already bound when extension reaches it; among those, prefer more C
    than G and a central position in the amplicon.
    """
    lo = forward.end + PROBE_PRIMER_GAP
    hi = reverse.start - PROBE_PRIMER_GAP
    if hi - lo < MIN_PROBE_LEN:
        return None

    highest_primer_tm = max(forward.tm, reverse.tm)
    tm_low = highest_primer_tm + PROBE_TM_ABOVE_PRIMER[0]
    tm_high = highest_primer_tm + PROBE_TM_ABOVE_PRIMER[1]
    tm_target = (tm_low + tm_high) / 2.0
    amplicon_centre = (forward.start + reverse.end) / 2.0

    best: Optional[Probe] = None
    best_key: Optional[Tuple[float, float, float]] = None
    for i in range(bisect.bisect_left(starts, lo), len(candidates)):
        probe = candidates[i]
        if probe.start > hi - MIN_PROBE_LEN:
            break
        if probe.end > hi:
            continue
        if not tm_low <= probe.tm <= tm_high:
            continue
        key = (-min(probe.c_minus_g, 0),
               abs(probe.tm - tm_target),
               abs((probe.start + probe.end) / 2.0 - amplicon_centre))
        if best_key is None or key < best_key:
            best_key, best = key, probe
    return best


# ============================================================================
# 7. ASSAY ASSEMBLY, SCORING AND RANKING
# ============================================================================

@dataclass
class Assay:
    target: str
    forward_primer: str
    forward_start: int          # 1-based, for reporting
    forward_tm: float
    forward_gc: float
    reverse_primer: str         # already reverse-complemented: order as synthesised
    reverse_start: int          # 1-based plus-strand start of its binding site
    reverse_tm: float
    reverse_gc: float
    probe: str
    probe_start: int            # 1-based
    probe_tm: float
    probe_gc: float
    amplicon_size: int
    cross_dimer: int
    score: float


def score_assay(forward: Oligo, reverse: Oligo, probe_tm: float,
                probe_gc: float, amplicon_size: int, cross_dimer: int) -> float:
    """Composite penalty score out of 100. Higher is better.

    Every term the design actually cares about is represented, including the
    probe. Scoring only the forward primer's GC makes most candidates tie.
    """
    optimal_primer_tm = (MIN_PRIMER_TM + MAX_PRIMER_TM) / 2
    probe_tm_target = max(forward.tm, reverse.tm) + sum(PROBE_TM_ABOVE_PRIMER) / 2

    score = 100.0
    score -= 2.0 * abs(forward.tm - optimal_primer_tm)
    score -= 2.0 * abs(reverse.tm - optimal_primer_tm)
    score -= 3.0 * abs(forward.tm - reverse.tm)          # matched Tm matters most
    score -= 0.2 * abs(forward.gc - 50.0)
    score -= 0.2 * abs(reverse.gc - 50.0)
    score -= 1.0 * abs(probe_tm - probe_tm_target)
    score -= 0.1 * abs(probe_gc - 55.0)
    score -= 0.05 * abs(amplicon_size - OPTIMAL_AMPLICON)
    score -= 1.5 * max(0, cross_dimer - MAX_3P_DIMER)
    return round(score, 2)


def design_assays(record: FastaRecord) -> List[Assay]:
    """Full design pass over one FASTA record."""
    sequence = record.sequence
    primers = find_candidate_primers(sequence)
    if not primers:
        return []

    probe_candidates, probe_starts = build_probe_candidates(sequence)
    if not probe_candidates:
        return []

    buckets = index_by_start(primers, len(sequence))
    assays: List[Assay] = []

    for forward in primers:
        # A reverse primer's binding site must start after the forward primer
        # ends, and the product must fall inside the amplicon window.
        first_start = max(forward.end, forward.start + MIN_AMPLICON - MAX_PRIMER_LEN)
        last_start = forward.start + MAX_AMPLICON - MIN_PRIMER_LEN
        for start in range(first_start, min(last_start, len(sequence)) + 1):
            for reverse in buckets[start]:
                amplicon_size = reverse.end - forward.start
                if not MIN_AMPLICON <= amplicon_size <= MAX_AMPLICON:
                    continue
                if abs(forward.tm - reverse.tm) > MAX_PRIMER_TM_DIFF:
                    continue

                reverse_oligo = reverse_complement(reverse.sequence)
                if three_prime_duplex_score(forward.sequence, reverse_oligo) > MAX_3P_DIMER:
                    continue
                if three_prime_duplex_score(reverse_oligo, forward.sequence) > MAX_3P_DIMER:
                    continue
                cross_dimer = max_duplex_score(forward.sequence, reverse_oligo)
                if cross_dimer > MAX_SELF_DIMER:
                    continue

                probe = select_probe(probe_candidates, probe_starts,
                                     forward, reverse)
                if probe is None:
                    continue

                assays.append(Assay(
                    target=record.identifier,
                    forward_primer=forward.sequence,
                    forward_start=forward.start + 1,
                    forward_tm=round(forward.tm, 2),
                    forward_gc=round(forward.gc, 1),
                    reverse_primer=reverse_oligo,
                    reverse_start=reverse.start + 1,
                    reverse_tm=round(reverse.tm, 2),
                    reverse_gc=round(reverse.gc, 1),
                    probe=probe.sequence,
                    probe_start=probe.start + 1,
                    probe_tm=round(probe.tm, 2),
                    probe_gc=round(probe.gc, 1),
                    amplicon_size=amplicon_size,
                    cross_dimer=cross_dimer,
                    score=score_assay(forward, reverse, probe.tm, probe.gc,
                                      amplicon_size, cross_dimer),
                ))
    return assays


def select_non_redundant(assays: List[Assay], limit: int = TOP_RESULTS,
                         separation: int = MIN_ASSAY_SEPARATION) -> List[Assay]:
    """Keep the best assay per region instead of ten shifts of one assay."""
    chosen: List[Assay] = []
    for assay in sorted(assays, key=lambda a: a.score, reverse=True):
        if all(assay.target != kept.target
               or abs(assay.forward_start - kept.forward_start) >= separation
               for kept in chosen):
            chosen.append(assay)
            if len(chosen) >= limit:
                break
    return chosen


# ============================================================================
# 8. OUTPUT
# ============================================================================

FIELDS = ["target", "forward_primer", "forward_start", "forward_tm", "forward_gc",
          "reverse_primer", "reverse_start", "reverse_tm", "reverse_gc",
          "probe", "probe_start", "probe_tm", "probe_gc",
          "amplicon_size", "cross_dimer", "score"]


def write_csv(assays: Iterable[Assay], path: str) -> None:
    rows = [asdict(a) for a in assays]
    try:
        import pandas as pd
        pd.DataFrame(rows, columns=FIELDS).to_csv(path, index=False)
    except ImportError:
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)


def write_oligo_fasta(assays: Sequence[Assay], path: str) -> None:
    """Write oligos as FASTA so they can be pasted into BLAST/Primer-BLAST."""
    with open(path, "w", encoding="utf-8") as handle:
        for i, a in enumerate(assays, start=1):
            handle.write(f">assay{i}_F {a.target} pos={a.forward_start}\n{a.forward_primer}\n")
            handle.write(f">assay{i}_R {a.target} pos={a.reverse_start}\n{a.reverse_primer}\n")
            handle.write(f">assay{i}_P {a.target} pos={a.probe_start}\n{a.probe}\n")


def print_table(assays: Sequence[Assay]) -> None:
    if not assays:
        print("\nNo assay met the criteria. Widen MIN/MAX_AMPLICON or the Tm "
              "window, or check the target for ambiguous bases.")
        return
    header = f"{'#':>2}  {'forward primer':<26}{'reverse primer':<26}{'probe':<32}" \
             f"{'Fwd Tm':>7}{'Rev Tm':>8}{'Prb Tm':>8}{'amp':>6}{'score':>8}"
    print("\n" + header)
    print("-" * len(header))
    for i, a in enumerate(assays, start=1):
        print(f"{i:>2}  {a.forward_primer:<26}{a.reverse_primer:<26}{a.probe:<32}"
              f"{a.forward_tm:>7.1f}{a.reverse_tm:>8.1f}{a.probe_tm:>8.1f}"
              f"{a.amplicon_size:>6}{a.score:>8.2f}")


# ============================================================================
# 9. ENTRY POINT
# ============================================================================

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="qPCR primer/probe design for pathogen detection.")
    parser.add_argument("--fasta", default=str(HERE / "pathogen.fasta"),
                        help="target FASTA (one or more records)")
    parser.add_argument("--out", default=str(HERE / "qpcr_assays.csv"),
                        help="output CSV path")
    parser.add_argument("--report-fasta", default=None,
                        help="also write the oligos as FASTA for BLAST screening")
    parser.add_argument("--top", type=int, default=TOP_RESULTS,
                        help="number of assays to report per run")
    parser.add_argument("--min-amplicon", type=int, default=MIN_AMPLICON)
    parser.add_argument("--max-amplicon", type=int, default=MAX_AMPLICON)
    parser.add_argument("--na", type=float, default=NA_CONC,
                        help="monovalent cation concentration, mol/L")
    parser.add_argument("--primer-conc", type=float, default=PRIMER_CONC,
                        help="primer concentration, mol/L")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global MIN_AMPLICON, MAX_AMPLICON, NA_CONC, PRIMER_CONC, PROBE_CONC
    args = parse_args(argv)
    MIN_AMPLICON, MAX_AMPLICON = args.min_amplicon, args.max_amplicon
    NA_CONC = args.na
    PRIMER_CONC = PROBE_CONC = args.primer_conc

    try:
        records = read_fasta(args.fasta)
    except FileNotFoundError:
        print(f"ERROR: no FASTA at {args.fasta}", file=sys.stderr)
        print(f"       Put pathogen.fasta in {HERE}, or pass --fasta <path>.",
              file=sys.stderr)
        return 1
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Loaded {len(records)} record(s) from {args.fasta}")
    for record in records:
        print("  " + summarise_record(record))

    all_assays: List[Assay] = []
    for record in records:
        found = design_assays(record)
        print(f"  {record.identifier}: {len(found)} valid assay(s)")
        all_assays.extend(found)

    selected = select_non_redundant(all_assays, limit=args.top)
    print_table(selected)

    if selected:
        write_csv(selected, args.out)
        print(f"\nWrote {len(selected)} assay(s) to {args.out}")
        if args.report_fasta:
            write_oligo_fasta(selected, args.report_fasta)
            print(f"Wrote oligos to {args.report_fasta}")
        print("REMINDER: specificity is not checked here. BLAST every oligo "
              "against the host and near-neighbour genomes before ordering.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:      # e.g. output piped into `head`
        sys.stderr.close()
        raise SystemExit(0)
