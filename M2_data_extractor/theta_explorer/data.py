"""Molecule payloads for the explorer: geometry, dipole and charges, via the package.

The explorer computes Sterimol and the dipole frame itself (ported, gated), but
the dipole vector and the partial charges come from the calculation, so they are
read here with DescriPyTor's own Molecule loader and passed through untouched.

    payload('mol.feather')  or  payload('mol.log')  ->  dict(name, xyz, dipole, charges)
"""
import contextlib
import io
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)                        # M2_data_extractor
for p in (PKG, os.path.dirname(PKG)):
    if p not in sys.path:
        sys.path.insert(0, p)


def _molecule(path):
    from data_extractor import Molecule
    path = os.path.abspath(path)
    if path.lower().endswith('.log'):               # Gaussian log -> feather, the package's own reader
        from feather_extractor import gauss_file_handler, save_to_feather
        tmp = os.path.join(tempfile.mkdtemp(prefix='theta_'), os.path.splitext(os.path.basename(path))[0])
        with contextlib.redirect_stdout(io.StringIO()):
            df, _ = gauss_file_handler(path)
            save_to_feather(df, tmp)
        path = tmp + '.feather'
    cwd = os.getcwd()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            return Molecule(path)
    finally:
        os.chdir(cwd)                                 # the loader changes directory


def payload(path, name=None):
    m = _molecule(path)
    x = m.xyz_df
    xyz = '%d\n%s\n' % (len(x), name or '') + '\n'.join(
        '%s %.8f %.8f %.8f' % (r.atom, r.x, r.y, r.z) for r in x.itertuples())
    dip = None
    d = getattr(m, 'gauss_dipole_df', None)
    if d is not None and len(d):
        dip = [float(v) for v in d.iloc[0].tolist()[:4]]           # dip_x, dip_y, dip_z, total
    charges = {}
    for t, df in (getattr(m, 'charge_dict', None) or {}).items():
        try:
            # the first n rows are the atoms (get_charge_df reads them the same way); NBO tables run on
            charges[t] = [float(v) for v in df.iloc[:len(x), 0].tolist()]
        except Exception:
            pass
    return dict(name=name or os.path.splitext(os.path.basename(path))[0], xyz=xyz, dipole=dip, charges=charges)


if __name__ == '__main__':
    # python data.py <file.feather|file.log> [fixture_name]
    #   with a fixture name, writes fixtures/<name>.data.json (dipole + charges) for the built-in examples
    import json
    p = payload(sys.argv[1])
    if len(sys.argv) > 2:
        name = sys.argv[2]
        out = os.path.join(HERE, 'fixtures', name + '.data.json')
        with open(out, 'w') as f:
            json.dump(dict(dipole=p['dipole'], charges=p['charges']), f)
        xyz = os.path.join(HERE, 'fixtures', name + '.xyz')      # the same atom order the charges are in
        if not os.path.exists(xyz):
            with open(xyz, 'w') as f:
                f.write(p['xyz'])
            print('wrote', xyz)
        print('wrote', out, '| dipole', p['dipole'], '| charges', sorted(p['charges']))
    else:
        print(json.dumps(dict(p, xyz=p['xyz'][:120] + '...'), indent=1)[:1200])
