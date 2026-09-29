"""SMILES -> numbered structures -> chained cluster stages (GOAT, UMA, xTB, ORCA, Gaussian).

    from M1_pre_calculations.pipeline import load, build_all, submit, status, retry, watch, fetch
    p = load("protocol.json")
    build_all(p); submit(p)

See docs/PIPELINE.md.
"""
from .build import BuildError, build_all, build_one
from .cluster import ClusterError, Remote, fetch, plan_lines, retry, status, submit, summarize, watch
from .protocol import Protocol, ProtocolError, load
from .stages import render_all, render_stage
