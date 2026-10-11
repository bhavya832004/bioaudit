# BioAudit methodology notes

These notes describe the similarity methods implemented in this repository. They are not a substitute for domain-specific homology criteria.

## Exact duplicates

Two records are exact duplicates when their normalized sequences are identical (uppercase, whitespace removed). Exact identity is reported separately from approximate similarity.

## `kmer_jaccard` (built-in)

- Metric name: `kmer_jaccard`.
- Definition: Jaccard index of the sets of contiguous k-mers from two sequences.
- This is **not** percent sequence identity, not an alignment, and not a homology test.
- The score depends on `k`, sequence length, and residue composition. Short sequences have few k-mers; compositionally biased sequences can inflate Jaccard overlap.
- Ambiguous or extended symbols (for example `X`) participate in k-mers as literal characters.
- Implementation: brute-force pairwise comparison. Default safety limit: 250,000 train×test pairs unless `--allow-quadratic` or a higher `--max-pairs` is set.
- An above-threshold score is a **risk indicator** that related sequences may appear in both partitions. It does not prove leakage or invalidate a published result.

## MMseqs2 sequence identity

- Metric name: `mmseqs_identity`.
- BioAudit calls `mmseqs easy-search` with `--min-seq-id`, `-c` coverage, and `--cov-mode`.
- Reported fields include identity, alignment length, and query/target coverage when present in tabular output.
- Sequence identity is **not** a universal homology criterion. Coverage mode and clustering policy matter and must be documented with any experiment.
- The `mmseqs` binary is optional. CI tests parse fixture TSV files and skip live searches when MMseqs2 is not installed.

## Interpretation

High train/test similarity can inflate apparent model performance. Detecting similar pairs supports a more careful split design; it does not by itself prove that a model or paper is invalid.
