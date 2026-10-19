"""The atom picker's pure helpers (M2_data_extractor/atom_picker.html), run in node: the page is one file, so each
function is cut out of it by name and called on its own."""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

PAGE = Path(__file__).resolve().parents[1] / "M2_data_extractor" / "atom_picker.html"
NODE = shutil.which("node")


def _function(src, name):
    """The source of `function name(...) {...}`, by brace count (these helpers hold no braces in strings)."""
    start = re.search(r"\bfunction %s\(" % re.escape(name), src).start()
    depth, i = 0, src.index("{", start)
    for j in range(i, len(src)):
        depth += {"{": 1, "}": -1}.get(src[j], 0)
        if depth == 0:
            return src[start:j + 1]
    raise ValueError(name)


def _run(names, body):
    if NODE is None:
        pytest.skip("node is not installed")
    src = PAGE.read_text(encoding="utf-8")
    script = "\n".join(_function(src, n) for n in names) + "\n" + body
    out = subprocess.run([NODE, "-e", script], capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr
    return json.loads(out.stdout)


def test_groups_are_checked_when_committed():
    got = _run(["groupProblem"], """
      const f = (key, size, label) => ({ key, size, label });
      console.log(JSON.stringify([
        groupProblem(f("sterimol", 2, "Sterimol  [origin, attached]"), [1, 4]),
        groupProblem(f("sterimol", 2, "Sterimol  [origin, attached]"), [1, 4, 5]),
        groupProblem(f("bond_angle", 0, "Bond angle/dihedral  [3 or 4 atoms]"), [8, 7, 4, 1]),
        groupProblem(f("bond_angle", 0, "Bond angle/dihedral  [3 or 4 atoms]"), [8, 7]),
      ]));""")
    assert got[0] == "" and got[2] == ""
    assert got[1] == "Sterimol takes 2 atoms per group, not 3."
    assert got[3] == "A bond angle takes 3 atoms and a dihedral 4, not 2."


def test_typed_atoms_read_as_groups():
    got = _run(["typedGroups"], """
      console.log(JSON.stringify([typedGroups("23,1,3 23,1,2,3"), typedGroups("1, 4   4 ,1"), typedGroups("1 4 7")]));""")
    assert got == [[["23", "1", "3"], ["23", "1", "2", "3"]], [["1", "4"], ["4", "1"]], [["1", "4", "7"]]]


def test_ring_picks_leave_the_page_one_ring_each():
    got = _run(["exportGroups"], """
      console.log(JSON.stringify([
        exportGroups({ key: "ring_vibration", size: 1 }, [[13], [20]]),
        exportGroups({ key: "charges", size: 1 }, [[1], [4]]),
        exportGroups({ key: "sterimol", size: 2 }, [[1, 4]]),
      ]));""")
    assert got == [[[13], [20]], [1, 4], [[1, 4]]]


def test_picks_made_on_another_molecule_are_named():
    got = _run(["picksThatChanged"], """
      console.log(JSON.stringify([
        picksThatChanged([4, 1, 9], ["C", "C", "O", "N"], ["C", "C", "O", "O"]),
        picksThatChanged([1, 2], [], ["C", "C"]),
      ]));""")
    assert got == [[4, 9], []]


def test_align_panel_reads_its_fields():
    got = _run(["alAtomRange", "alAnchors", "alEls", "alStalled"], """
      console.log(JSON.stringify([
        alAtomRange("1-4 9, 12 x"), alAnchors("m_Br: 9, 9\\nbad line\\no_F:3 4"), alEls("3\\nc\\nC 0 0 0\\nH 1 0 0\\nBr 0 1 0"),
        alStalled({ rmsd: 6.0, rmsd_kabsch: 0.49 }), alStalled({ rmsd: 0.5, rmsd_kabsch: 0.49 }), alStalled(null),
      ]));""")
    assert got == [[1, 2, 3, 4, 9, 12], {"m_Br": [[9, 9]], "o_F": [[3, 4]]}, ["C", "H", "Br"], True, False, False]


def test_a_swap_exchanges_two_numbers_everywhere():
    got = _run(["alSwapAtoms", "alEls"], """
      const r = { old_to_new: { 1: 2, 2: 1, 3: 3 }, xyz: "3\\nx\\nO 0 0 0\\nC 1 0 0\\nH 2 0 0", stage: ["a", "b", "c"], differs: [] };
      alSwapAtoms(r, 2, 3, ["C", "O", "H"]);
      console.log(JSON.stringify([r.old_to_new, r.xyz.split("\\n").slice(2), r.stage, r.differs, r.edited]));""")
    assert got == [{"1": 3, "2": 1, "3": 2}, ["O 0 0 0", "H 2 0 0", "C 1 0 0"], ["a", "c", "b"], [1, 2, 3], True]


def test_a_field_the_files_cannot_give_says_why():
    got = _run(["fieldMissing"], """
      const NEEDS = { dipole: "dipole", charges: "charges", ring_vibration: "vibrations" };
      let INPUTS = { shared: { dipole: false, vibrations: true, charges: [] } };
      const a = [fieldMissing("sterimol"), fieldMissing("dipole"), fieldMissing("charges"), fieldMissing("ring_vibration")];
      INPUTS = null;
      a.push(fieldMissing("dipole"));
      console.log(JSON.stringify(a));""")
    assert got == ["", "the files in this folder have no dipole", "the files in this folder have no charges", "", ""]


def test_morfeus_sterimol_pairs_stay_with_morfeus():
    got = _run(["pickedSterimolPairs"], """
      const state = { sterimol: [[1, 4]], mf_sterimol: [[2, 7]] };
      console.log(JSON.stringify(pickedSterimolPairs()));""")
    assert got == [[1, 4]]
