"""
Reference and alternate DNA sequences of a variant in a fixed-length window, built as in the variant scoring
of the training data, which used kipoiseq 0.7.1 the way the Enformer usage notebook does:

    interval  = Interval(chrom, start, end).resize(length)
    anchor    = interval.center() - interval.start          # relative; kipoiseq clamps it into the window
    reference = VariantSeqExtractor(...).extract(interval, [], anchor=anchor)
    alternate = VariantSeqExtractor(...).extract(interval, [variant], anchor=anchor)

For one variant the alternate is the reference window with REF replaced by ALT; deletions extend and
insertions truncate the window on the side away from the anchor, so its length is unchanged. REF is not
checked against the genome (1_prepare_variants.py orients every variant to hg38 first). This module
reproduces that construction with pyfaidx alone, so the scoring environments do not need kipoiseq;
featurization/README.md describes the test against kipoiseq.

    genome = Genome("hg38.fa")
    start, end = window(pos - 1, pos - 1, 393216)        # Enformer: centred on the variant
    ref_seq, alt_seq = sequences(genome, "chr11", start, end, pos, "A", "G")
"""
import numpy as np
from pyfaidx import Fasta, Sequence


class Genome:
    """Reference genome; positions outside a chromosome read as N (the Enformer notebook's FastaStringExtractor)."""

    def __init__(self, fasta_file):
        self.fasta = Fasta(fasta_file)
        self.sizes = {name: len(record) for name, record in self.fasta.items()}

    def fetch(self, chrom, start, end):
        """Uppercase sequence of the 0-based half-open interval [start, end)."""
        size = self.sizes[chrom]
        s, e = max(start, 0), min(end, size)
        seq = str(self.fasta.get_seq(chrom, s + 1, e).seq).upper() if e > s else ""
        return "N" * max(-start, 0) + seq + "N" * max(end - size, 0)


def window(start, end, width):
    """(start, end) of kipoiseq's Interval(chrom, start, end).resize(width) on the + strand."""
    if end - start == width:
        return start, end
    center = (start + end) // 2 + (start + end) % 2
    return center - width // 2 - width % 2, center + width // 2


def sequences(genome, chrom, start, end, pos, ref, alt):
    """Reference and alternate sequence of the window [start, end) for the variant chrom:pos:ref:alt
    (pos 1-based); the steps of kipoiseq's VariantSeqExtractor.extract with fixed_len=True."""
    reference = genome.fetch(chrom, start, end)
    anchor = (start + end) // 2 + (start + end) % 2 - start      # what the notebook passes (relative)
    anchor = max(min(anchor, end), start)                        # kipoiseq clamps it into the window
    s = pos - 1
    pairs = [(Sequence(name=chrom, seq=ref, start=s, end=s + len(ref)),
              Sequence(name=chrom, seq=alt, start=s, end=s + len(alt)))]
    split = []                                                   # 1. split a variant spanning the anchor
    for r, a in pairs:
        if r.start < anchor < r.end:
            mid = anchor - r.start
            split += [(r[:mid], a[:mid]), (r[mid:], a[mid:])]
        else:
            split.append((r, a))
    up = sorted([p for p in split if p[0].start >= anchor], key=lambda p: p[0].start)
    down = sorted([p for p in split if p[0].start < anchor], key=lambda p: p[0].start, reverse=True)
    istart, iend = start, end                                    # 2. extend the window for deletions
    for r, a in up:
        if len(a) < len(r):
            iend += len(r) - len(a)
    for r, a in down:
        if len(a) < len(r):
            istart -= len(r) - len(a)
    seq = genome.fetch(chrom, istart, iend)

    def ref_piece(x, y):
        i = x - istart
        return seq[i: i + max(0, y - x)]

    down_parts, prev = [], anchor                                # 3. build outwards from the anchor
    for r, a in down:
        if r.end <= istart:
            break
        down_parts += [ref_piece(r.end, prev), a.seq]
        prev = r.start
    down_parts.append(ref_piece(istart, prev))
    up_parts, prev = [], anchor
    for r, a in up:
        if r.start >= iend:
            break
        up_parts += [ref_piece(prev, r.start), a.seq]
        prev = r.end
    up_parts.append(ref_piece(prev, iend))
    down_str = "".join(reversed(down_parts))
    up_str = "".join(up_parts)
    down_len, up_len = anchor - start, end - anchor              # 4. cut to the window length
    down_str = down_str[-down_len:] if down_len else ""
    up_str = up_str[:up_len] if up_len else ""
    return reference, down_str + up_str


def one_hot(seq):
    """(length, 4) float32 one-hot over A, C, G, T; N = 0.25 each, other characters 0
    (kipoiseq.transforms.functional.one_hot_dna, used by the Enformer notebook)."""
    table = np.zeros((256, 4), dtype=np.float32)
    table[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.eye(4, dtype=np.float32)
    table[ord("N")] = 0.25
    return table[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
