"""
Command line driver for cylindrical ESP sampling, called by esp_grid_job.sh.

Split into a write step and a read step because xtb runs between them. The
geometry is rotated into the Sterimol frame BEFORE xtb sees it, and the rotated
geometry is what gets submitted: the potential is a scalar field, so rotating
the molecule does not change it, and having xtb and the grid share one frame
removes any chance of the grid being expressed in different axes than the
structure it is meant to sample.
"""

import argparse
import pickle

import pandas as pd

from M2_data_extractor.extractor_utils.esp_grid_utils import (
    angular_moments,
    build_cylindrical_grid,
    read_xtb_esp,
    write_esp_coord,
)
from M2_data_extractor.extractor_utils.sterimol_utils import (
    preform_coordination_transformation,
)
from utils.help_functions import get_df_from_file


def write_step(args) -> None:
    xyz_df = get_df_from_file(args.xyz)
    xyz_df[['x', 'y', 'z']] = xyz_df[['x', 'y', 'z']].astype(float)

    base_atoms = [int(i) for i in args.base_atoms.split(',')] if args.base_atoms else None
    framed = preform_coordination_transformation(xyz_df, indices=base_atoms,
                                                 origin=args.origin)

    grid = build_cylindrical_grid(framed, n_phi=args.n_phi, n_axial=args.n_axial,
                                  probe_offset=args.probe_offset, radii=args.radii)
    n_points = write_esp_coord(grid, args.coord)

    # The frame the grid lives in is the frame xtb must be given.
    with open(args.framed_xyz, 'w') as handle:
        handle.write(f'{len(framed)}\n')
        handle.write('rotated into the Sterimol frame for ESP sampling\n')
        for _, row in framed.iterrows():
            handle.write(f'{row["atom"]:<3s} {row["x"]:>15.8f} '
                         f'{row["y"]:>15.8f} {row["z"]:>15.8f}\n')

    with open(args.grid, 'wb') as handle:
        pickle.dump(grid, handle)

    print(f'wrote {n_points} grid points to {args.coord} '
          f'({args.n_axial} levels x {args.n_phi} azimuths)')


def read_step(args) -> None:
    with open(args.grid, 'rb') as handle:
        grid = pickle.load(handle)

    # Raises if the point count or the coordinates disagree, which is how a
    # stale esp_coord surfaces as an error instead of as quiet wrong numbers.
    sampled = read_xtb_esp(grid, args.esp)
    sampled.to_csv(args.out, index=False)

    moments = angular_moments(sampled, beta=args.beta)
    moments.to_csv(args.moments, index=False)

    # |Z| near zero means the field has no angular structure to align, so it is
    # worth seeing before any descriptor built on the phases is believed.
    summary = moments[['Z_S_1', 'Z_S_2', 'Z_E_1', 'Z_E_2']].abs().mean()
    print(f'sampled {len(sampled)} points; mean |Z|: '
          + ', '.join(f'{k}={v:.3f}' for k, v in summary.items()))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='step', required=True)

    write = sub.add_parser('write', help='build the grid and write esp_coord')
    write.add_argument('--xyz', required=True)
    write.add_argument('--grid', default='grid.pkl')
    write.add_argument('--coord', default='esp_coord')
    write.add_argument('--framed-xyz', default='framed.xyz')
    write.add_argument('--base-atoms', default=None,
                       help='comma separated indices defining the Sterimol frame')
    write.add_argument('--origin', default=None)
    write.add_argument('--n-phi', type=int, default=36)
    write.add_argument('--n-axial', type=int, default=20)
    write.add_argument('--probe-offset', type=float, default=0.5)
    write.add_argument('--radii', default='bondi', choices=['bondi', 'CPK'])
    write.set_defaults(func=write_step)

    read = sub.add_parser('read', help='attach xtb_esp.dat to the grid')
    read.add_argument('--grid', default='grid.pkl')
    read.add_argument('--esp', default='xtb_esp.dat')
    read.add_argument('--out', required=True)
    read.add_argument('--moments', required=True)
    read.add_argument('--beta', type=float, default=1.0)
    read.set_defaults(func=read_step)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
