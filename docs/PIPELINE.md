# Cluster pipeline: SMILES to finished structures

`descripytor pipeline` takes a CSV of SMILES and does five things:

1. builds starting geometries with a fixed atom numbering;
2. writes one SGE job script per stage (ORCA GOAT, UMA, GFN2-xTB, ORCA DFT, Gaussian);
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
descripytor pipeline retry   protocol.json            # resubmit transient failures from their first failed stage
descripytor pipeline reset   protocol.json --only m001,m006 --from-stage 2   # forget markers, to run those stages again
descripytor pipeline fetch   protocol.json            # last stage's outputs -> <name>_pipeline/fetched/
descripytor pipeline adopt   protocol.json --stage 1 --sources sources.csv   # outside results as a finished stage
descripytor pipeline submit  protocol.json --from-stage 2 --only-adopted     # continue from the next stage
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

Every molecule's stages are separate SGE jobs, each held on that molecule's previous stage
(`qsub -hold_jid`), so nothing polls and one slow molecule never holds up another. The molecule
reaches the script as `-v PIPE_TASK=<line in ids.txt>`. The BGU cluster runs N1GE 6.0u8, which
has no `-terse`, no per-task `-hold_jid_ad`, and rejects `qsub -t`, so the chain is built by
`submit_chain.sh` from a plan file. `jobs.json` records every job id. Every job does the same
checks:

- it reads the previous stage's `out/<id>.xyz`;
- it refuses the input if the element order differs from the manifest;
- it runs in node scratch, falling back to `$TMPDIR`;
- it writes `status/<id>.<stage>`: `running`, `done`, or `failed: <reason>`.

A job whose molecule is already `done` at that stage exits at once.

| kind | does | output | options (defaults) |
|---|---|---|---|
| `goat` | ORCA GOAT conformer search | lowest conformer; the ensemble in `work/` | `keywords` (`! GOAT XTB`), `maxcore` |
| `uma` | FAIRChem UMA reoptimisation (LBFGS) | reoptimised xyz; the comment line says whether it converged | `fmax` (0.05), `steps` (300), `model` (`uma-s-1p1`), `task` (`omol`) |
| `xtb` | GFN2-xTB optimisation (with the Hessian, `--ohess`), then a single point | xyz, `.q` charges, `.wbo`, `.props` (dipole, HOMO, LUMO), `.vibspectrum` in `work/` | `level` (`--gfn 2`), `opt_level` (xtb's `crude` … `vtight` … `extreme`; default xtb's `normal`), `hess` (true), `imag_tol` (20), `imag_retry` (1; 0 = never displace along a mode), `metal_contacts` (`allow`; `constrain`: when the free optimisation puts the metal on the ligand, e.g. Ni(0) side-on to an arene, re-optimise from the input with every M–D–X angle fixed at its input value, so the metal stays σ-bound on the donor axis; labelled in `checks/`) |
| `orca` | any ORCA job; fails on imaginary frequencies | optimised xyz; `.out`, `.property.txt` and `.hess` in `work/` | `keywords` (`! r2SCAN-3c Opt Freq`), `maxcore` |
| `gaussian` | Gaussian 16; must be the last stage | `out/<id>.log`, ready for `descripytor logs_to_feather` | `route` (**required**), `mem_gb` (32), `tail` (text after the coordinates) |

All stages also take `cores` and `resources` (a list of extra `#$ -l` lines, e.g.
`["h_vmem=4G"]`). The queue is `cluster.queue`, the same for every stage.

## Checks: what makes a stage fail

A stage is `done` only when every one of these checks passes. When one fails, the output is
moved to `out/rejected/` and the next stage cannot start: each stage requires the molecule's
previous stage to be marked `done`, not just an output file to exist.

| check | stages | fails when |
|---|---|---|
| program errors | GOAT, ORCA | "aborting the run", "error termination", no "ORCA TERMINATED NORMALLY"; ORCA also fails on an optimisation or SCF that did not converge |
| | xTB | "abnormal termination", "failed to converge", `#ERROR`, "convergence criteria cannot be satisfied", no "normal termination of xtb", an optimisation that did not converge |
| | UMA | a Python error, or the step limit reached without convergence (`allow_unconverged: true` turns this into a warning) |
| | Gaussian | "Error termination", or no "Normal termination" at the end of the log |
| imaginary frequencies | xTB (`hess`, on by default: `--ohess`), ORCA and Gaussian frequency jobs | any mode below −`imag_tol` (default 20 cm⁻¹). Smaller ones are logged as warnings. xTB first restarts `imag_retry` times (default 1) from its geometry displaced along the mode (`xtbhess.xyz`), the usual cure for a saddle point. |
| structure | every stage that writes an xyz, and `adopt` | a non-metal bond broken or formed relative to the build (package bond rule: 1.15 × Pyykkö radii); two atoms closer than 0.6 × their radii sum; a donor or ancillary more than 1.25 × (r_M + r_X) from the metal; a chelate bite outside 65–105° |

Warnings don't stop a molecule. They go to `checks/<id>.<stage>`, show as `done*` in `status`, and
`status` lists them. Examples: a new metal contact such as an arene or agostic interaction, a GOAT
worker that did not converge, an imaginary mode removed by a restart.

**Retrying:** `retry` resubmits transient failures only, such as a killed job or a node problem.
A rejected structure, imaginary frequencies or an unconverged UMA run would repeat with the same
input, so those need a change first. `retry --all` resubmits them anyway.

## Adopting finished outside results

Structures computed before a pipeline existed can join it as a completed stage. The typical case
is GOAT minima from a hand-written job. `sources.csv` has three columns:

- `id`: the pipeline id (`m001`, …);
- `remote_source`: the finished file on the cluster;
- `local_start`: the outside run's input xyz, relative to the CSV.

**How it works:**

1. The outside run may use another atom order. The permutation is found by matching the
   pipeline's build to `local_start` atom by atom on coordinates. This works only when both were
   built identically; anything else is refused, never guessed.
2. On the cluster, `adopt.sh` reorders each source into `<stage>/out/<id>.xyz`, checks the element
   order, and marks the stage done.
3. Sources that are not there yet are reported as `missing` and left alone.

**Continuing:** `submit --from-stage k+1 --only-adopted` submits the adopted molecules from the next
stage. `jobs.json` is merged, and a molecule already submitted at that stage is refused. Re-run
`adopt` and `submit` as the outside jobs finish.

## What runs where

- **On the cluster:** only bash on the login node (`status.sh`, `qsub`); UMA and xTB run on the
  compute nodes inside the jobs.
- **Transfers:** base64 over ssh, because a login banner corrupts scp. The tcsh login shell
  swallows ssh's standard input and rejects any command-line word over about 8,000 characters, so
  an upload goes as 4,000-character parts in numbered files, joined remotely and checked by md5.
- **`watch`** is a loop on your machine; if the machine sleeps, run `status` later instead.

## Remote layout

```
<root>/<name>/
  protocol.json  manifest.json  ids.txt  jobs.json  plan.txt  status.sh  submit_chain.sh
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
