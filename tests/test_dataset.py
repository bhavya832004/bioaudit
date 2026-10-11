import json
from urllib.error import URLError

import pytest

from bioaudit.dataset import (
    download_uniprot_fasta,
    generate_toy_classification_records,
    write_toy_classification_csv,
)
from bioaudit.io import read_records
from bioaudit.provenance import sha256_file


class _DummyResponse:
    def __init__(self, payload: bytes):
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self):
        return self._payload


def test_toy_classification_is_offline_and_labelled(tmp_path):
    out = tmp_path / "toy.csv"
    path, manifest = write_toy_classification_csv(out)
    records = read_records(path, label_column="label")
    assert len(records) == 21
    labels = {r["label"] for r in records}
    assert labels == {"nucleus", "cytoplasm", "membrane"}
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["source"] == "toy_classification"
    assert payload["sha256"] == sha256_file(path)
    assert payload["record_count"] == 21
    families = {r["family"] for r in generate_toy_classification_records()}
    assert len(families) == 9


def test_uniprot_download_uses_injected_opener(tmp_path):
    fasta = ">sp|P00001|X Test\nACDEFGHIKLMNPQRSTVWY\n>sp|P00002|Y Test\nGAVLIPFYWSTCMNQDEKRH\n"

    def opener(request, timeout=90):
        assert timeout == 5
        assert "uniprot.org" in request.full_url
        return _DummyResponse(fasta.encode("utf-8"))

    out = tmp_path / "uni.fasta"
    path, manifest = download_uniprot_fasta(2, out, timeout=5, opener=opener)
    records = read_records(path)
    assert len(records) == 2
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["record_count"] == 2
    assert payload["suitable_supervised_benchmark"] is False
    assert payload["requested_limit"] == 2


def test_uniprot_rejects_non_fasta(tmp_path):
    with pytest.raises(RuntimeError, match="no FASTA"):
        download_uniprot_fasta(
            10,
            tmp_path / "bad.fasta",
            opener=lambda *a, **k: _DummyResponse(b'{"error": "nope"}'),
        )


def test_uniprot_retries_then_fails(tmp_path):
    def opener(request, timeout=90):
        raise URLError("network down")

    with pytest.raises(RuntimeError, match="UniProt download failed"):
        download_uniprot_fasta(10, tmp_path / "x.fasta", opener=opener, retries=1)


def test_uniprot_limit_bounds(tmp_path):
    with pytest.raises(ValueError, match="limit"):
        download_uniprot_fasta(0, tmp_path / "x.fasta")
