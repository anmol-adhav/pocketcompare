"""Unit tests for fetch module — no network required."""
import pytest
from pathlib import Path
from pocketcompare.fetch import load_local


def test_load_local_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_local(tmp_path / "nonexistent.pdb")


def test_load_local_exists(tmp_path):
    f = tmp_path / "test.pdb"
    f.write_text("ATOM  ...")
    result = load_local(f)
    assert result == f.resolve()
