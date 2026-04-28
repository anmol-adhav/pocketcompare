"""Fetch structures from RCSB PDB, AlphaFold DB, or load local files."""
from pathlib import Path
import tempfile
import requests

RCSB_URL = "https://files.rcsb.org/download/{code}.pdb"
AFDB_URL  = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"


def fetch_pdb(code: str, out_dir: str | Path | None = None) -> Path:
    """Download a PDB structure from RCSB by 4-letter code."""
    code = code.upper()
    out_dir = Path(out_dir) if out_dir else Path(tempfile.mkdtemp())
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{code}.pdb"
    resp = requests.get(RCSB_URL.format(code=code), timeout=60)
    resp.raise_for_status()
    path.write_bytes(resp.content)
    return path


def fetch_alphafold(uniprot_id: str, out_dir: str | Path | None = None) -> Path:
    """Download an AlphaFold v4 model by UniProt ID."""
    out_dir = Path(out_dir) if out_dir else Path(tempfile.mkdtemp())
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"AF-{uniprot_id}.pdb"
    resp = requests.get(AFDB_URL.format(uniprot=uniprot_id), timeout=60)
    resp.raise_for_status()
    path.write_bytes(resp.content)
    return path


def load_local(path: str | Path) -> Path:
    """Resolve and validate a local PDB file path."""
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"File not found: {p}")
    return p
