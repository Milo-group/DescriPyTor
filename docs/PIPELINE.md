# Cluster pipeline: SMILES to finished structures

`descripytor pipeline` takes a CSV of SMILES and does five things:

1. builds starting geometries with a fixed atom numbering;
2. writes one SGE array script per stage (ORCA GOAT, UMA, GFN2-xTB, ORCA DFT, Gaussian);
3. uploads the run;
4. chains the stages so each molecule moves on the moment its previous step ends;
5. reports status, retries failures and fetches the results.

It is driven by one protocol file, which is also the record of how the structures were made.
Code: `M1_pre_calculations/pipeline/`.

```bash
descripytor pipeline build   protocol.json            # SMILES -> numbered xyz, manifest
descripytor pipeline submit  protocol.json --dry-run  # render the scripts, show the qsub chain
descripytor pipeline submit  protocol.json            # upload, submit every stage at once
descripytor pipeline status  protocol.json            # one row per molecule, one column per stage
descripytor pipeline watch   protocol.json            # report every 10 min until finished
descripytor pipeline retry   protocol.json            # resubmit each molecule from its first failed stage
descripytor pipeline fetch   protocol.json            # last stage's outputs -> <name>_pipeline/fetched/
```

## The protocol

```json
{
  "name": "phosphines_ni",
  "molecules": "phosphines.csv",
  "smiles_column": "smiles",
  "name_column": "ligand",
  "charge": 0,
  "multiplicity": 1,
  "build": {"type": "metal_mono", "metal": "Ni", "donor_smarts": "[PX3]", "metal_distance": 2.2},
  "stages": [
    {"kind": "goat", "cores": 8},
    {"kind": "uma", "cores": 8},
    {"kind": "xtb", "cores": 4},
    {"kind": "gaussian", "cores": 16, "mem_gb": 32, "route": "#p M062X/def2TZVP opt freq pop=(nbo,hirshfeld)"}
  ],
  "cluster": {"host": "...", "key": "~/.ssh/id_ed25519", "root": "/gpfs0/.../pipelines", "queue": "fairshare.q",
              "scratch": "/local/gaus", "orca": "...", "mpi": "...", "g16root": "...",
              "uma_env": "...", "hf_token_file": "..."}
}
```

A complete example for the BGU cluster ships as
`M1_pre_calculations/pipeline/examples/protocol_bgu.json`. Relative paths are read from the
protocol file's folder.

## Building: the numbering is decided here

| `build.type` | atom order | options |
|---|---|---|
| `organic` | with `core_smarts`: the core atoms first, in SMARTS order, so a descriptor atom has the same number in every molecule | `core_smarts` |
| `metal_mono` | metal, donor, then the ligand | `metal`, `donor_smarts` (default `[PX3]`, must pick one atom), `metal_distance` |
| `metal_chelate` | metal, ancillaries, donor 1, donor 2, then the ligand: the order `MetalComplex` reads | `metal`, `ancillary` (`H` or `Cl`), `n_ancillary` (2: square planar, trans to the donors; 1: trigonal), `donor_smarts` with atom maps `:1` and `:2` (otherwise found from the graph), `metal_distance` (a number, or a dict per donor element) |

All types take `n_conformers` (default 10) and `seed` (default 0).

**How a molecule is built:**

1. RDKit embeds `n_conformers` conformers (ETKDG), relaxes them with MMFF (UFF where MMFF lacks
   parameters), and tries them in energy order.
2. The metal goes on the donor lone pair.
3. A structure is accepted when no atom other than the donors sits within 2.2 Å (heavy) or 1.8 Å
   (H) of the metal, and, for a chelate, when the bite angle is between 65° and 105°.

**Chelates:**

- The two donors are held 2.5–3.0 Å apart while embedding, because a free ligand such as
  bipyridine otherwise comes out anti.
- The ligand is then relaxed around the fixed metal, so a methyl or tBu hydrogen turns out of
  the pocket. The metal placement and ancillary positions follow CS3's `build_general.py`.
- A monodentate ligand gets the same relaxation only when no conformer leaves room for the metal.

**Outputs, in `<name>_pipeline/`:**

- `build/<id>.xyz`;
- `elements/<id>.elements`: the element order every later stage is checked against;
- `ids.txt`;
- `manifest.json`: id, name, SMILES, formula, and the metal, donor, core and ancillary indices,
  1-based.

Ids are `m001`, `m002`, … unless the CSV has an `id_column`.

## Stages

Every stage is one SGE array over the molecules. Task *i* of stage *k* waits only for task *i* of
stage *k−1* (`qsub -hold_jid_ad`), so nothing polls. Every task does the same checks:

- it reads the previous stage's `out/<id>.xyz`;
- it refuses the input if the element order differs from the manifest;
- it runs in node scratch, falling back to `$TMPDIR`;
- it writes `status/<id>.<stage>`: `running`, `done`, or `failed: <reason>`.

A task that is already `done` exits at once.

| kind | does | output | options (defaults) |
|---|---|---|---|
| `goat` | ORCA GOAT conformer search | lowest conformer; the ensemble in `work/` | `keywords` (`! GOAT XTB`), `maxcore` |
| `uma` | FAIRChem UMA reoptimisation (LBFGS) | reoptimised xyz; the comment line says whether it converged | `fmax` (0.05), `steps` (300), `model` (`uma-s-1p1`), `task` (`omol`) |
| `xtb` | GFN2-xTB optimisation, then a single point | xyz, `.q` charges, `.wbo`, `.props` (dipole, HOMO, LUMO) | `level` (`--gfn 2`) |
| `orca` | any ORCA job; fails on imaginary frequencies | optimised xyz; `.out`, `.property.txt` and `.hess` in `work/` | `keywords` (`! r2SCAN-3c Opt Freq`), `maxcore` |
| `gaussian` | Gaussian 16; must be the last stage | `out/<id>.log`, ready for `descripytor logs_to_feather` | `route` (**required**), `mem_gb` (32), `tail` (text after the coordinates) |

All stages also take `cores` and `resources` (a list of extra `#$ -l` lines, e.g.
`["h_vmem=4G"]`). The queue is `cluster.queue`, the same for every stage.

## What runs where

- **On the cluster:** only bash on the login node (`status.sh`, `qsub`); UMA and xTB run on the
  compute nodes inside the jobs.
- **Transfers:** base64 over ssh, because a login banner corrupts scp. Uploads go in 20,000-character
  chunks and are checked by md5.
- **`watch`** is a loop on your machine; if the machine sleeps, run `status` later instead.

## Remote layout

```
<root>/<name>/
  protocol.json  manifest.json  ids.txt  jobs.json  status.sh
  build/<id>.xyz   elements/<id>.elements   status/<id>.<stage>
  s1_goat/  run.sh  logs/  work/<id>.out  out/<id>.xyz
  s2_uma/   run.sh  uma_one.py  logs/  work/  out/
  ...
```

## After fetching

For metal complexes, read `fetched/s*_xtb/out/` with `MetalComplex`, pairing each xyz with its
`.q`, `.wbo` and `.props` in `XtbSinglePoint`. For the Gaussian route, run
`descripytor logs_to_feather` on `fetched/s*_gaussian/out/`, then `Molecules`. The numbering in
`manifest.json` is the numbering to use in the input JSON. Then follow
[EXTRACTION_PROTOCOL.md](EXTRACTION_PROTOCOL.md).
