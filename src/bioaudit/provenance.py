"""Dataset and run provenance helpers."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def environment_metadata() -> dict[str, Any]:
    return {
        "bioaudit_version": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "executable": sys.executable,
    }


def write_manifest(path: str | Path, payload: dict[str, Any]) -> Path:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    document = dict(payload)
    document.setdefault("generated_at_utc", utc_now_iso())
    document.setdefault("environment", environment_metadata())
    dest.write_text(json.dumps(document, indent=2, default=str), encoding="utf-8")
    return dest
