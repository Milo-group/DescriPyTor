# DescriPyTor — install and GUI (one page)

Two ways. **Docker is the easy one** (chemistry software is already inside the image). Use pip only if you already live in conda.

**Visual version** (diagrams on GitHub, or open in a browser and print to PDF): [docs/visual-guide.md](docs/visual-guide.md) · [docs/visual-guide.html](docs/visual-guide.html)

---

## A. Docker (recommended)

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) and wait until it says the engine is running.
2. Get the code and start the apps:

```bash
git clone https://github.com/Milo-group/DescriPyTor.git
cd DescriPyTor
mkdir work
docker compose up --build
```

First build takes several minutes. Leave that terminal open.

3. In a browser:

| What | Open |
|---|---|
| 3D atom picker + modeling | http://localhost:7432/visual |
| Streamlit extraction app | http://localhost:8503 |

4. Drop molecule files (`.feather`, `.xyz`, CSVs) into the `work` folder on your computer. Inside the GUI, use paths like `/work/yourfile.feather`. Results written to `/work` show up in that same folder.

Stop: `Ctrl+C` in the terminal, or `docker compose down`.

---

## B. pip (no Docker)

Needs a conda env with Python 3.10 or 3.11.

```bash
conda create -n descripytor python=3.10
conda activate descripytor
pip install "descripytor[gui,webapp]"
```

**Desktop GUI** (older Tk window):

```bash
descripytor gui
```

**Browser GUI** (same as Docker, from a clone):

```bash
git clone https://github.com/Milo-group/DescriPyTor.git
cd DescriPyTor
python M2_data_extractor/gui_server.py
```

Then open http://localhost:7432/visual

---

## Use it in five minutes

Example molecules ship with the package (26 substituted benzenes).

**Python**

```python
from M2_data_extractor import Molecules
from descripytor.examples import feather_example_dir, input_example_json

mols = Molecules(str(feather_example_dir()), threshold=1.82)
print(mols.success_molecules)
```

**Command line** (from a clone, so the JSON path is local):

```bash
descripytor extractor -i Getting_started_with_examples/input_example.json -o features -f Getting_started_with_examples/feather_example
```

Or, after `pip install`, copy the JSON out of the package:

```python
from pathlib import Path
from descripytor.examples import input_example_json, feather_example_dir
print(input_example_json())
print(feather_example_dir())
```

**In the picker:** load a `.feather` or `.xyz` → click atoms → extract → optional model on the same page. Atom numbers are **1-based** (Gaussian style).

Gaussian `.log` files → feathers first:

```bash
descripytor logs_to_feather
```

---

## If something fails

| Symptom | Fix |
|---|---|
| `docker compose` cannot connect | Start Docker Desktop and wait until it is idle |
| Port 7432 or 8503 already in use | Close the other app, or set `DESCRIPYTOR_GUI_PORT` / `DESCRIPYTOR_HOST_PORT` |
| `descripytor gui` missing customtkinter | `pip install "descripytor[gui]"` |
| RDKit / Morfeus import errors on pip | Use Docker instead |

Longer docs: [README.md](README.md).
