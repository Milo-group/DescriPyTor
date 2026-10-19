# Conformaty test case 03: ring numbering

Copied from Conformaty (`test_cases/03_ring_numbering`, MIT licence, https://github.com/Milo-group/Conformaty).
`stage_order_ref.xyz` is cyclohexane (chair); `stage_order_query.xyz` is a methyl-substituted oxane, so the query has a
ring O where the reference has a carbon, and one more carbon.

More than one numbering is defensible, and the order of the look-up steps decides which one you get:

| steps | ring O | methyl C |
|---|---|---|
| default (complete rings, extend the common part, any element by position) | 6 | 8 |
| complete rings only | 6 | 8 |
| extend the common part only | 19 (numbered last) | 6 |

Checked against Conformaty's own `pair` (default order: same heavy-atom numbering, 2026-10-19). Each run takes under a
minute on a laptop. In the visual extractor, *Align & renumber* > *Try Conformaty's example* runs it on a copy of this
folder; reorder the steps under *More options* and run it again to see the numbering change.
