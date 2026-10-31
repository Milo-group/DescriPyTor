"""Forwards to feature_explorer/data.py: this folder is the explorer's name before it was renamed,
kept so scripts that put it on sys.path and `import data` still run."""
import importlib.util, os, sys
_spec = importlib.util.spec_from_file_location(__name__, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "feature_explorer", "data.py"))
_mod = importlib.util.module_from_spec(_spec); sys.modules[__name__] = _mod; _spec.loader.exec_module(_mod)
