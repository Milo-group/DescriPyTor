"""Where things are. Inside the package nothing has to be configured; the environment can override each entry.

  descripytor      the package this folder is in: its feature explorer is the molecule view, its readers read the structures
  engine           this folder (engine.py and core_fixed.py: the sampler set-up every run uses)
  bassa_python     the Python running the server; new tables need `bassa_reg` in it (`descripytor bassa_explorer --bassa-python`)
  reader_python    BASSA_EXPLORER_READER: the Python that has pyarrow and the package's readers (default: the one that started the command)
  cases            BASSA_EXPLORER_CASES: the folder where kept cases live (default: ./bassa_cases); never inside the package
  paper, cs2_feathers   BASSA_EXPLORER_PAPER, BASSA_EXPLORER_CS2_FEATHERS: only for rebuilding the prepared cases of the DescriPyTor paper
  chrome           CHROME: only for test_page.py
"""
import os, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
E = os.environ.get
S = dict(bassa_python=sys.executable, reader_python=E("BASSA_EXPLORER_READER") or sys.executable, engine=str(HERE), descripytor=str(HERE.parents[1]),
         paper=E("BASSA_EXPLORER_PAPER", ""), cs2_feathers=E("BASSA_EXPLORER_CS2_FEATHERS", ""),
         chrome=E("CHROME") or "C:/Program Files/Google/Chrome/Application/chrome.exe")
CASES_DIR = Path(E("BASSA_EXPLORER_CASES") or Path.cwd() / "bassa_cases").resolve()
