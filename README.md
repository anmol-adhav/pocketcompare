# pocketcompare

A minimal tool to compare protein binding pockets and support engineering decisions.

## Two modes

| Mode | Question answered |
|---|---|
| `homolog` | Given two proteins, should I engineer for **ACTIVITY** or **FLUX**? |
| `mutant` | Did my mutation **expand**, **contract**, or leave the pocket **unchanged**? |

## Install

```bash
git clone https://github.com/yourusername/pocketcompare
cd pocketcompare
pip install -e .
```

## Homolog mode

```bash
pocketcompare homolog \
  --pdb1 2ZCP --pdb2 5IYS \
  --pairs "A:172,A:173;A:176,A:177;A:225,A:237"
```

## Mutant mode

```bash
pocketcompare mutant \
  --wt 2ZCP \
  --mut W368F_model.pdb \
  --pocket "A:172;A:176;A:225;A:228" \
  --wt-label WT --mut-label W368F
```

## Find pocket residues

```bash
python examples/find_pocket_residues.py 2ZCP
```

## Run tests

```bash
# fast unit tests only (no network)
pytest -m "not integration"

# all tests including network calls
pytest
```

## Decision guide

### Homolog mode
| Conservation ratio | Decision |
|---|---|
| < 0.50 | **ACTIVITY** — mutate pocket residues |
| 0.50–0.80 | **ACTIVITY first**, then test flux |
| > 0.80 | **FLUX** — pathway-level changes more impactful |

### Mutant mode
| Effect | What it means | Next step |
|---|---|---|
| EXPANDED | Less steric clash, better substrate/product access | Whole-cell flux assay |
| CONTRACTED | Tighter fit, may improve Km | In vitro kinetics first |
| UNCHANGED | Acts via electrostatics/dynamics | MD simulation or B-factor analysis |
