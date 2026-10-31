"""Forwards to M2_data_extractor/feature_packages.py. The module moved into the package so that an installed DescriPyTor has it
(the visual extractor's Extract CSV runs it); this file keeps `import feature_packages` and `python feature_packages.py ...` working from here."""
import importlib.util, os, sys
_spec = importlib.util.spec_from_file_location(__name__, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "M2_data_extractor", "feature_packages.py"))
_mod = importlib.util.module_from_spec(_spec); sys.modules[__name__] = _mod; _spec.loader.exec_module(_mod)
