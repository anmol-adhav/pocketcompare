# 1. unzip and install
pip install -e ".[dev]"

# 2. run fast unit tests first (no network needed)
pytest -m "not integration" -v

# 3. run integration tests (fetches real PDB structures)
pytest -m integration -v

# 4. try the examples
python examples/find_pocket_residues.py 2ZCP
python examples/homolog_example.py
python examples/mutant_example.py
