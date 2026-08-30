# DescriPyTor — install and GUI (one page)

Start with **conda + pip**. Use Docker only if that fails, or you need xTB / extra engines.

**Visual version** (diagrams on GitHub, or open in a browser and print to PDF): [docs/visual-guide.md](docs/visual-guide.md) · [docs/visual-guide.html](docs/visual-guide.html)

---

## A. conda + pip (recommended)

Needs conda and Python 3.10 or 3.11.

```bash
git clone https://github.com/Milo-group/DescriPyTor.git
cd DescriPyTor
conda create -n descripytor python=3.10
conda activate descripytor
pip install -e .
python M2_data_extractor/gui_server.py
```

Then open http://localhost:7432/visual

If RDKit or igraph fail on pip:

```bash
conda install -c conda-forge rdkit python-igraph
```

Then run the GUI server again.

**Optional desktop window** (older Tk app): `pip install -e ".[gui]"` then `descripytor gui`.

---

## Use it in five minutes

Example molecules ship with the package (26 substituted benzenes).

**In the picker:** load a `.feather` or `.xyz` → click atoms → extract → optional model on the same page. Atom numbers are **1-based** (Gaussian style). Caffeine loads as a demo until you choose a file.

**Python**

```python
from M2_data_extractor import Molecules
from descripytor.examples import feather_example_dir

mols = Molecules(str(feather_example_dir()), threshold=1.82)
print(mols.success_molecules)
```

**Command line** (from the clone):

```bash
descripytor extractor -i Getting_started_with_examples/input_example.json -o features -f Getting_started_with_examples/feather_example
```

Gaussian `.log` files → feathers first:

```bash
descripytor logs_to_feather
```

---

## B. Docker (if conda fails, or you need the full toolkit)

Chemistry software is already in the image (RDKit, Morfeus, xTB, Streamlit). Heavier first step: Docker Desktop, then a source build.

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) and wait until it says the engine is running.
2. From the clone:

```bash
cd DescriPyTor
mkdir work
docker compose up --build
```

First build takes several minutes. Leave that terminal open.

3. In a browser:

| What | Open |
|---|---|
| 3D atom picker | http://localhost:7432/visual |
| Feature extraction GUI | http://localhost:7432/ |
| Streamlit extraction app | http://localhost:8503 |

4. Drop molecule files into `work/` on your computer. Inside the GUI, use paths like `/work/yourfile.feather`. Results written to `/work` show up in that same folder.

Stop: `Ctrl+C` in the terminal, or `docker compose down`.

---

## If something fails

| Symptom | Fix |
|---|---|
| RDKit / igraph / Morfeus import errors | `conda install -c conda-forge rdkit python-igraph`, then Docker if that still fails |
| Port 7432 already in use | Close the other app, or set `GUI_PORT` |
| `docker compose` cannot connect | Start Docker Desktop and wait until it is idle |
| Port 8503 already in use | Close the other app, or set `DESCRIPYTOR_HOST_PORT` |

Longer docs: [README.md](README.md).
