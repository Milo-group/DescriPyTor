# Small test set (8 alcohols + 1 tiny Ni complex)

Self-contained XYZ geometries for extraction and modeling tests. No Gaussian
logs, feathers, or xTB dumps are required.

## Files

| Path | What it is |
|---|---|
| `xyz/*.xyz` | One conformer per molecule |
| `targets.csv` | Name, SMILES, boiling point (°C) |
| `modeling_table.csv` | Extracted Sterimol / bonds / %Vbur + `bp_c` (written by `extract.py`) |
| `classification_toy.csv` | 12-row subset of the logistic example (class labels) |

## Atom numbering (every alcohol)

1. C<sub>α</sub> (the carbon bound to OH)
2. O
3. H of OH
4. C<sub>β</sub> when present

Sterimol of the alkyl group uses origin **2 (O)** and direction **1 (C)**.

Geometries are tetrahedral Z-matrices (`build_xyz.py`), not DFT minima. Bond
lengths and Sterimol *trends* (methanol < n-butanol in L; tert-butanol larger
B5 than methanol) are what the tests check, not published crystal structures.

## Rebuild

From the repository root:

```powershell
python tests/data/small_set/build_xyz.py
python tests/data/small_set/extract.py
```

`extract.py` uses `MetalComplex` Sterimol / %Vbur (the CS3 kernel), not morfeus.
Topology columns (`tf_*`) are added when RDKit is installed.
