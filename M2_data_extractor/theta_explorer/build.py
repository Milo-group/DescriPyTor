r"""Assemble the Sterimol theta explorer into one self-contained HTML page.

sterimol.js (a line-for-line port of get_sterimol_df) and scene.js (the renderer)
are inlined into template.html with the fixture molecules, so the page runs with
no server and no network. gui_server.py serves assemble() at /theta; running
this file writes a standalone theta_explorer.html next to it.

    node gate.js && python build.py
    python build.py --inject ..\..\Getting_started_with_examples\descriptor_extraction_toolkit\atom_picker.html
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# Paper presets. Atom numbers are 1-based, as in the package; the columns they
# reproduce are named so a preset can be checked against its table. view='axis'
# opens a preset in the theta figure's own view of its active axis.
_FLAGS = ('axis', 'L', 'b1', 'plane', 'loc', 'b5', 'theta', 'drop', 'inset')
_SHOW = lambda **k: {**dict.fromkeys(_FLAGS, False), 'axis': True, 'drop': True, **k}            # noqa: E731
_ALL = dict.fromkeys(_FLAGS, True)
_DIP = dict(on=True, origin='23-28', y='23', plane='1')
_NONE = dict(axes=[], active=-1, dip=dict(on=False, origin='', y='', plane='', show=dict(frame=True, mu=True, u=True, v=False, w=False)),
             sites=[], map='', fade='', fadeOff=False, legend=True)
_L = lambda **k: dict(_NONE, **k)                                                                 # noqa: E731
# the case-study figure panels: the CYLview look of the paper's other renders, labels as symbols (values are in the plots)
PANEL_SCALE = 38.6   # drawing units per Angstrom: one scale across the three case-study panels
# the published look: every hydrogen, small balls and hard outlines, so the constructions
# read through the substituents; black bonds keep the element colours for the atoms alone
_PANEL = dict(look='paper', hydrogens='all', atomLabels='off', values=False, labelScale=1.35,
              preset='glossy', palette='paper', bondColor='black', outline=1.6,
              atomScale=0.54, bondScale=0.44, hiddenLines=True, shadow=False, fog=0)


# Figure 2: the three constructions the paper adds, all on one substrate so they can be read
# against each other -- the dipole in a frame the molecule defines, the place along the axis
# where B5 is read, and the angle between v and the B1 plane. Quiet ink, no leader lines.
_F2 = dict(look='paper', ink='colour', preset='glossy', palette='paper', bondColor='black',
           hydrogens='all', atomLabels='off', values=True, featureLabels=True, leaders=False,
           labelScale=2.0, outline=1.4, atomScale=0.5, bondScale=0.4, hiddenLines=True,
           shadow=False, fog=0, locMark='tick')
_F2_SCALE = 34.0
# the opening figure: the same three constructions with the numbers taken out -- symbols only,
# no value key, no plane caption. What each feature IS, before any of them has a value.
_F1 = dict(_F2, values=False, planeNote=False, labelScale=2.1)
# the substituent ladder under the opening figure: one axis per group, drawn alone at one scale
_GRP = dict(_F2, labelScale=2.4, hydrogens='all', planeNote=False)
_GRP_SCALE = 41.0


def _group(name, a, b):
    return dict(name='Group: %s' % name, view='axis', perA=_GRP_SCALE, style=_GRP,
                layers=_L(axes=[dict(a=a, b=b, color='#151A21',
                                     show=_SHOW(L=True, b1=True, plane=True, loc=True, b5=True, theta=True))],
                          active=0, legend=False, fadeOff=True, fadeHide=True))


# The three panels of the feature figure, in the paper's own idiom: one construction each,
# everything the construction does not touch dropped, quiet type, no leader lines.
_FIG = dict(look='paper', ink='colour', preset='glossy', palette='paper', bondColor='black',
            hydrogens='all', atomLabels='off', values=True, featureLabels=True, leaders=False,
            labelScale=1.9, outline=1.4, atomScale=0.5, bondScale=0.4, hiddenLines=True,
            shadow=False, fog=0.15, locMark='tick', planeNote=False)


def _theta(a, b):
    return dict(name='theta construction on %d -> %d: B1 and its plane, v and its shadow, loc B5, B5, end-on view' % (a, b),
                view='axis', layers=_L(axes=[dict(a=a, b=b, color='#1596A6', show=_ALL)], active=0, fadeOff=True))


FIXTURES = [
    dict(file='p-Otf.xyz', name='p-OTf', sub='Case Study 1 substrate', a=1, b=23, presets=[
        _theta(1, 23),
        dict(name='Case Study 1: mu_u, q_CM5 on C1, q_NBO on C3 (dipole_y_{23..28}-23-1, cm5_atom_1, nbo_atom_3)',
             layers=_L(dip=dict(_DIP, show=dict(frame=False, mu=True, u=True, v=False, w=False)),
                       sites=[dict(i=1, type='cm5', color='#C9352B'), dict(i=3, type='nbo', color='#E58A1F')])),
        dict(name='Fig. 2B: the dipole frame u, v, w and all three projections',
             layers=_L(dip=dict(_DIP, show=dict(frame=True, mu=True, u=True, v=True, w=True)))),
        dict(name='Fig. 2A: Sterimol C1 -> C23 with B1, loc B5, B5 and theta',
             layers=_L(axes=[dict(a=1, b=23, color='#1596A6', show=_SHOW(b1=True, loc=True, b5=True, theta=True))], active=0)),
        dict(name='CM5 charge map', layers=_L(map='cm5')),
        dict(name='Figure 2a: the dipole in the frame the molecule defines (mu, u, mu_u)',
             view=dict(right=[27, 25], up=[23, 26], perA=_F2_SCALE), style=_F2,
             layers=_L(dip=dict(_DIP, show=dict(frame=False, mu=True, u=True, v=False, w=False),
                                colors=dict(mu='#98A0AB', u='#7A3E9D')), legend=False)),
        dict(name='Figure 1a: mu and its projection, symbols only',
             view=dict(right=[27, 25], up=[23, 26], perA=_F2_SCALE), style=_F1,
             layers=_L(dip=dict(_DIP, show=dict(frame=False, mu=True, u=True, v=False, w=False),
                                colors=dict(mu='#98A0AB', u='#7A3E9D')), legend=False)),
        dict(name='Figure 1b: L, B1, B5 and loc B5, symbols only',
             view=dict(right=[1, 23], up=[23, 26], perA=_F2_SCALE), style=_F1,
             layers=_L(axes=[dict(a=1, b=23, color='#151A21',
                                  show=_SHOW(L=True, b1=True, loc=True, b5=True))], active=0, legend=False, fadeOff=True)),
        dict(name='Figure 2b: L, B1, B5 and the loc B5 they are read at (C1 -> C23)',
             view=dict(right=[1, 23], up=[23, 26], perA=_F2_SCALE), style=_F2,
             layers=_L(axes=[dict(a=1, b=23, color='#151A21',
                                  show=_SHOW(L=True, b1=True, loc=True, b5=True))], active=0, legend=False, fadeOff=True)),
        dict(name='Figure 2c: theta, between v and the B1 plane',
             view='axis', style=_F2,
             layers=_L(axes=[dict(a=1, b=23, color='#151A21',
                                  show=_SHOW(b1=True, plane=True, b5=True, theta=True))], active=0, legend=False, fadeOff=True)),
    ]),
    dict(file='m-CN.xyz', name='m-CN', sub='Case Study 1 substrate (the published panel)', a=1, b=23, presets=[
        dict(name='Case Study 1 panel: mu_u along ring centroid -> ipso, q_CM5 on C1, q_NBO on C3 (dipole_y_{23..28}-23-1, cm5_atom_1, nbo_atom_3)',
             view=dict(right=[27, 25], up=[23, 26], perA=PANEL_SCALE), style=_PANEL,
             layers=_L(dip=dict(_DIP, show=dict(frame=False, mu=True, u=True, v=False, w=False), colors=dict(mu='#98A0AB', u='#7A3E9D')),
                       sites=[dict(i=1, type='cm5', color='#C9352B'), dict(i=3, type='nbo', color='#E58A1F')], legend=False)),
        dict(name='Figure 2a: the dipole in the frame the molecule defines (mu, u, mu_u)',
             view=dict(right=[27, 25], up=[23, 26], perA=_F2_SCALE), style=_F2,
             layers=_L(dip=dict(_DIP, show=dict(frame=False, mu=True, u=True, v=False, w=False),
                                colors=dict(mu='#98A0AB', u='#7A3E9D')), legend=False)),
        dict(name='Feature a: the dipole in the ring frame, against the raw vector',
             view=dict(right=[27, 25], up=[23, 26], perA=33.0), style=dict(_FIG, values=False),
             layers=_L(dip=dict(_DIP, rawAxes=True, show=dict(frame=True, mu=True, u=True, v=False, w=False),
                                colors=dict(mu='#7A3E9D', u='#1D64C4', v='#C9352B', w='#2E9E5B')),
                       fade='2-22', fadeHide=True, legend=True)),
        dict(name='Feature b: loc B5, where along L the widest point sits',
             view=dict(right=[1, 23], up=[23, 26], perA=44.0), style=_FIG,
             layers=_L(axes=[dict(a=1, b=23, color='#151A21',
                                  show=_SHOW(L=True, b1=True, loc=True, b5=True))], active=0,
                       legend=False, fadeOff=True, fadeHide=True)),
        _group('methyl', 1, 2),
        _group('benzyl', 8, 9),
        _group('phenyl', 9, 12),
dict(name='Figure 1a: mu and its projection, symbols only',
             view=dict(right=[27, 25], up=[23, 26], perA=_F2_SCALE), style=_F1,
             layers=_L(dip=dict(_DIP, show=dict(frame=False, mu=True, u=True, v=False, w=False),
                                colors=dict(mu='#98A0AB', u='#7A3E9D')), legend=False)),
        dict(name='Figure 1b: L, B1, B5 and loc B5, symbols only',
             view=dict(right=[1, 23], up=[23, 26], perA=_F2_SCALE), style=_F1,
             layers=_L(axes=[dict(a=1, b=23, color='#151A21',
                                  show=_SHOW(L=True, b1=True, loc=True, b5=True))], active=0, legend=False, fadeOff=True)),
        dict(name='Figure 2b: L, B1, B5 and the loc B5 they are read at (C1 -> C23)',
             view=dict(right=[1, 23], up=[23, 26], perA=_F2_SCALE), style=_F2,
             layers=_L(axes=[dict(a=1, b=23, color='#151A21',
                                  show=_SHOW(L=True, b1=True, loc=True, b5=True))], active=0, legend=False, fadeOff=True)),
        dict(name='Figure 2c: theta, between v and the B1 plane',
             view='axis', style=_F2,
             layers=_L(axes=[dict(a=1, b=23, color='#151A21',
                                  show=_SHOW(b1=True, plane=True, b5=True, theta=True))], active=0, legend=False, fadeOff=True)),
        _theta(1, 23),
    ]),
    dict(file='FL_lig_13_str_target.xyz', name='FL_lig_13', sub='Case Study 2 ligand', a=6, b=7, presets=[
        dict(name='Case Study 2 panel: loc B5,1 and B1 on R1 (6->4), loc B5,2 on R2 (6->7) (loc_B5_6-4, B1_6-4, loc_B5_6-7)',
             view=dict(right=[4, 6], up=[15, 5], perA=PANEL_SCALE), style=_PANEL,
             layers=_L(axes=[dict(a=6, b=4, color='#E07A1F', tag='1', show=_SHOW(b1=True, loc=True)),
                             dict(a=6, b=7, color='#C9352B', tag='2', show=_SHOW(loc=True))], active=0, legend=False)),
        _theta(6, 7),
    ]),
    dict(file='045_lig.xyz', name='045_lig', sub='Case Study 3 CuCl complex (tBu-BOX)', a=8, b=12,
         note='Case Study 3 reads these arms with MetalComplex: sub_B1_sym 2.751 &Aring; and sub_theta_sym 38.7&deg; in the '
              'model&rsquo;s table, against 2.754 &Aring; and 38.3&deg; here.', presets=[
        dict(name='Feature c: theta, how far the bulk leans off the B1 plane',
             view='axis', perA=76.0, style=dict(_FIG, values=False),
             layers=_L(axes=[dict(a=8, b=12, color='#151A21',
                                  show=_SHOW(b1=True, plane=True, b5=True, theta=True))], active=0,
                       legend=True, fadeOff=True, fadeHide=True)),
        _group('tert-butyl', 8, 12),
        dict(name='Figure 2c: theta, between v and the B1 plane (the case study 3 arm)',
             view='axis', style=_F2,
             layers=_L(axes=[dict(a=8, b=12, color='#151A21',
                                  show=_SHOW(b1=True, plane=True, b5=True, theta=True))], active=0, legend=False, fadeOff=True)),
        dict(name='Case Study 3 panel: on each arm the axis from C*, B1, and theta between v and the B1 plane (sub_B1_sym, sub_theta_sym)',
             view=dict(right=[3, 4], up=[1, [3, 4]], yaw=-40, perA=PANEL_SCALE), style=_PANEL,
             layers=_L(axes=[dict(a=8, b=12, color='#1596A6', show=_SHOW(b1=True, theta=True, drop=False)),
                             dict(a=19, b=21, color='#1596A6', show=_SHOW(b1=True, theta=True, drop=False))], active=0, legend=False)),
        _theta(8, 12),
    ]),
    dict(file='041_lig.xyz', name='041_lig', sub='Case Study 3 Cu complex', a=7, b=6,
         note='Case Study 3 descriptors come from MetalComplex, which measures from the same stereocentre but '
              'excludes it and starts its scan from a different frame: 1.8923 &Aring; and 27.85&deg; there, '
              'against the platform values shown here.', presets=[
        _group('isopropyl', 7, 6),
        dict(name='Figure 1c: theta, symbols only',
             view='axis', perA=82.0, style=_F1,
             layers=_L(axes=[dict(a=7, b=6, color='#151A21',
                                  show=_SHOW(b1=True, plane=True, b5=True, theta=True))], active=0, legend=False,
                       fadeOff=True, fadeHide=True)),
        dict(name='Figure 2c: theta, between v and the B1 plane (the case study 3 arm)',
             view='axis', perA=82.0, style=_F2,
             layers=_L(axes=[dict(a=7, b=6, color='#151A21',
                                  show=_SHOW(b1=True, plane=True, b5=True, theta=True))], active=0, legend=False,
                       fadeOff=True, fadeHide=True)),
        dict(name='Case Study 3: the isopropyl arm (7 -> 6), B1, the B5 direction and theta',
             layers=_L(axes=[dict(a=7, b=6, color='#1596A6', show=_SHOW(b1=True, b5=True, theta=True))], active=0)),
        _theta(7, 6),
    ]),
]


def _read(*parts):
    with open(os.path.join(HERE, *parts), encoding='utf-8') as f:
        return f.read()


# The same panel on three more members of each set. Case study 1's substrates and case study 2's
# ligands are built in one atom order, so the panel's numbers carry across them unchanged (checked:
# the element sequence is identical up to the substituent). Case study 3's complexes are not, so the
# two arms are found by substructure instead -- [#8]1[#6][#6]([!#1])[#7]~[#6]1, whose match gives the
# ring carbon and the group it carries. On 045 that returns 8->12 and 19->21, the published pair.
_CS1_PANEL = dict(view=dict(right=[27, 25], up=[23, 26], perA=PANEL_SCALE), style=_PANEL,
                  layers=_L(dip=dict(_DIP, show=dict(frame=False, mu=True, u=True, v=False, w=False), colors=dict(mu='#98A0AB', u='#7A3E9D')),
                            sites=[dict(i=1, type='cm5', color='#C9352B'), dict(i=3, type='nbo', color='#E58A1F')], legend=False))
_CS2_PANEL = dict(view=dict(right=[4, 6], up=[15, 5], perA=PANEL_SCALE), style=_PANEL,
                  layers=_L(axes=[dict(a=6, b=4, color='#E07A1F', tag='1', show=_SHOW(b1=True, loc=True)),
                                  dict(a=6, b=7, color='#C9352B', tag='2', show=_SHOW(loc=True))], active=0, legend=False))


def _cs3_panel(arm1, arm2):
    return dict(view=dict(right=[3, 4], up=[1, [3, 4]], yaw=-40, perA=PANEL_SCALE), style=_PANEL,
                layers=_L(axes=[dict(a=arm1[0], b=arm1[1], color='#1596A6', show=_SHOW(b1=True, theta=True, drop=False)),
                                dict(a=arm2[0], b=arm2[1], color='#1596A6', show=_SHOW(b1=True, theta=True, drop=False))],
                          active=0, legend=False))


# Figure 1: one parameter, several positions. The cone and the sphere are drawn from the
# position given here; Sterimol uses the axis the fixture declares.
_POS_PANEL = dict(style=dict(_PANEL, values=True, labelScale=1.3),
                  layers=_L(axes=[], active=-1))


# P-C runs 1.83-1.87 A, past the package's flat 1.82 A cutoff, so the phosphines declare their own
PHOS_THR = 1.95
MORE = [
    ('026_box_CuCl.xyz', 'Cu-BOX', 'bis(oxazoline) complex', 1, 3, _POS_PANEL),
    # a %Vbur ladder on one scaffold, 48 to 77 % at the metal, over two metals: something to
    # push the cone and the sphere against. All xTB single structures from the CS3 set.
    ('Ni_BOX_072.xyz', 'Ni-BOX 072', 'bis(oxazoline), smallest', 1, 4, _POS_PANEL),
    ('Ni_BOX_078.xyz', 'Ni-BOX 078', 'bis(oxazoline)', 1, 4, _POS_PANEL),
    ('Cu_PyOx_041.xyz', 'Cu-PyOx 041', 'pyridine-oxazoline', 1, 3, _POS_PANEL),
    ('Cu_BOX_018.xyz', 'Cu-BOX 018', 'bis(oxazoline)', 1, 3, _POS_PANEL),
    ('Cu_BOX_047.xyz', 'Cu-BOX 047', 'iPr-bis(oxazoline)', 1, 3, _POS_PANEL),
    ('Cu_BOX_020.xyz', 'Cu-BOX 020', 'bis(oxazoline)', 1, 3, _POS_PANEL),
    ('Cu_BOX_008.xyz', 'Cu-BOX 008', 'bis(oxazoline), bulkiest', 1, 3, _POS_PANEL),
    ('045_box_CuCl.xyz', 'Cu-tBuBOX', 'bis(oxazoline) complex', 1, 3, _POS_PANEL),
    # the Hammett benzoic acids, staged from the study-case feathers, so each carries its
    # dipole and three charge sets as well: the axis is ring -> substituent
    ('p-H_BA.xyz', 'benzoic acid', 'benzoic acid', 1, 7, _POS_PANEL),
    ('p-F_BA.xyz', 'p-F benzoic acid', 'benzoic acid', 4, 11, _POS_PANEL),
    ('p-Cl_BA.xyz', 'p-Cl benzoic acid', 'benzoic acid', 4, 11, _POS_PANEL),
    ('p-OMe_BA.xyz', 'p-OMe benzoic acid', 'benzoic acid', 4, 12, _POS_PANEL),
    ('p-COMe_BA.xyz', 'p-COMe benzoic acid', 'benzoic acid', 4, 12, _POS_PANEL),
    ('p-CO2Me_BA.xyz', 'p-CO2Me benzoic acid', 'benzoic acid', 1, 7, _POS_PANEL),
    ('p-t-Bu_BA.xyz', 'p-tBu benzoic acid', 'benzoic acid', 4, 12, _POS_PANEL),
    ('m-OMe_BA.xyz', 'm-OMe benzoic acid', 'benzoic acid', 3, 12, _POS_PANEL),
    ('m-CN_BA.xyz', 'm-CN benzoic acid', 'benzoic acid', 3, 12, _POS_PANEL),
    ('m-SO2Me_BA.xyz', 'm-SO2Me benzoic acid', 'benzoic acid', 3, 12, _POS_PANEL),
    ('o-Me_BA.xyz', 'o-Me benzoic acid', 'benzoic acid', 2, 11, _POS_PANEL),
    ('naphthyl_BA.xyz', '2-naphthoic acid', 'benzoic acid', 5, 11, _POS_PANEL),
    ('009_phos_CuCl.xyz', 'Cu-PHOX', 'phosphine complex', 1, 3, _POS_PANEL),
    ('009_phos.xyz', 'PHOX (free)', 'phosphine ligand', 15, 19, _POS_PANEL),
    ('055_phos.xyz', 'PCy-PHOX (free)', 'phosphine ligand', 10, 11, _POS_PANEL),
    ('p-NO2.xyz', 'p-NO2', 'Case Study 1 substrate', 1, 23, _CS1_PANEL),
    ('m-OMe.xyz', 'm-OMe', 'Case Study 1 substrate', 1, 23, _CS1_PANEL),
    ('p-CF3.xyz', 'p-CF3', 'Case Study 1 substrate', 1, 23, _CS1_PANEL),
    ('FL_lig_1.xyz', 'FL_lig_1', 'Case Study 2 ligand', 6, 7, _CS2_PANEL),
    ('FL_lig_7.xyz', 'FL_lig_7', 'Case Study 2 ligand', 6, 7, _CS2_PANEL),
    ('FL_lig_12.xyz', 'FL_lig_12', 'Case Study 2 ligand', 6, 7, _CS2_PANEL),
    ('026_lig.xyz', '026_lig', 'Case Study 3 CuCl complex (Me-BOX)', 6, 5, _cs3_panel((6, 5), (17, 19))),
    ('032_lig.xyz', '032_lig', 'Case Study 3 CuCl complex (iPr-BOX)', 7, 6, _cs3_panel((7, 6), (17, 19))),
    ('046_lig.xyz', '046_lig', 'Case Study 3 CuCl complex (aryl-BOX)', 8, 12, _cs3_panel((8, 12), (22, 24))),
]
FIXTURES += [dict(file=f, name=n, sub=s, a=a, b=b,
                  **({'thr': PHOS_THR} if 'phos' in f else {}),
                  presets=[dict(panel, name='%s panel: the published construction on %s' % (s.split(' ligand')[0].split(' substrate')[0].split(' CuCl')[0], n)),
                           _theta(a, b)])
             for f, n, s, a, b, panel in MORE]


def assemble():
    mols = []
    for m in FIXTURES:
        m = dict(m, key=m['file'], xyz=_read('fixtures', m['file']))
        data = os.path.join(HERE, 'fixtures', os.path.splitext(m['file'])[0] + '.data.json')
        if os.path.isfile(data):                       # dipole + charges from the calculation (data.py)
            with open(data) as f:
                m.update(json.load(f))
        mols.append(m)
    refs = {}
    for c in json.loads(_read('refs.json')):
        refs.setdefault(c['file'], {})['%d-%d' % (c['a'], c['b'])] = c['ref']
    scene = _read('scene.js')
    scenecss = scene.split('const SCENE_CSS = `', 1)[1].split('`;', 1)[0]
    return (_read('template.html')
            .replace('/*STERIMOL*/', _read('sterimol.js'))
            .replace('/*SCENE*/', scene)
            .replace('/*SMARTS*/', _read('smarts.js'))
            .replace('/*PROMOL*/', _read('promol.js'))
            .replace('/*VBUR*/', _read('vbur.js'))
            .replace('/*SCENECSS*/', scenecss)
            .replace('/*MOLECULES*/', json.dumps(mols))
            .replace('/*REFS*/', json.dumps(refs)))


BEGIN, END = '<!-- THETA-EXPLORER:BEGIN', '<!-- THETA-EXPLORER:END -->'


def inject(path):
    """Refresh the explorer embedded in an atom-picker page, between its markers.

    The page is stored as a JS string with every '<' escaped, so nothing in it
    (</script>, <!--) can end the host page's script early.
    """
    with open(path, encoding='utf-8') as f:
        html = f.read()
    i, j = html.index(BEGIN), html.index(END)
    payload = json.dumps(assemble()).replace('<', r'\u003c')
    block = (BEGIN + ' generated by M2_data_extractor/theta_explorer/build.py --inject; do not edit -->\n'
             '<script>window.THETA_EXPLORER_HTML = ' + payload + ';</script>\n')
    with open(path, 'w', encoding='utf-8', newline='') as f:
        f.write(html[:i] + block + html[j:])
    return len(payload)


if __name__ == '__main__':
    import sys
    if sys.argv[1:2] == ['--inject']:
        for page in sys.argv[2:]:
            print('%s  +%.0f KB explorer' % (page, inject(page) / 1024))
        sys.exit(0)
    html = assemble()
    out = os.path.join(HERE, 'theta_explorer.html')
    with open(out, 'w', encoding='utf-8') as f:
        f.write(html)
    print('%s  %.0f KB' % (out, len(html.encode()) / 1024))
