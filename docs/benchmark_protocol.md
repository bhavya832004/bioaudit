# Benchmark protocol

This protocol is independent of any single dataset.

1. Prepare data with `scripts/prepare_dataset.py` and keep the manifest.
2. Validate sequences and labels (`bioaudit audit`).
3. Cluster sequences (default: k-mer Jaccard connected components; MMseqs2 clustering is optional later).
4. Build a random stratified split and a homology-aware split. Clusters must not cross partitions in the homology-aware split.
5. Fit character k-mer logistic regression on training sequences only.
6. Evaluate on the test partition with macro-F1, balanced accuracy, per-class scores, and a confusion matrix.
7. Repeat for predeclared seeds. Do not tune on test.
8. Report the random minus homology-aware macro-F1 difference in percentage points, with caveats.

The default config (`configs/experiment.yaml`) uses the synthetic toy fixture. Those numbers are pipeline evidence only, not a biological result.
