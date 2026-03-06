#!/usr/bin/env python3
"""
Compute the structural distance between two extended XYZ configurations.

The distance is the RMSD (root-mean-square displacement) between matching
atoms.  Atoms are matched by their 'id' field when present, otherwise by
line order.  Periodic boundary conditions (minimum image convention) are
applied using the lattice from file A.

Outputs:
  - RMSD (Å)
  - Max per-atom displacement and which atom
  - RMS per component (dx, dy, dz)
  - Energy difference (if both files carry an energy)

For trajectory files, use --frame-a / --frame-b to select frames.

Usage:
  ./xyz_distance.py A.xyz B.xyz
  ./xyz_distance.py traj.xyz B.xyz --frame-a 5
  ./xyz_distance.py A.xyz B.xyz --no-pbc
  ./xyz_distance.py A.xyz B.xyz --verbose       # print per-atom displacements
"""

import argparse
import math
import re
import sys
from collections import OrderedDict


# ── XYZ parser (shared with xyz2lammps.py) ──────────────────────────────────

def _parse_properties(prop_str):
    parts = prop_str.strip().split(':')
    result = []
    for i in range(0, len(parts), 3):
        result.append((parts[i].lower(), parts[i+1].upper(), int(parts[i+2])))
    return result


def _parse_comment(line):
    info = {}
    m = re.search(r'[Ll]attice\s*=\s*"([^"]*)"', line)
    if m:
        info['lattice'] = list(map(float, m.group(1).split()))
    m = re.search(r'[Pp]roperties\s*=\s*(\S+)', line)
    if m:
        info['properties'] = _parse_properties(m.group(1))
    m = re.search(r'[Ee]nergy\s*[=:]\s*(-?[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)', line)
    if m:
        info['energy'] = float(m.group(1))
    return info


def read_frame(lines, frame_idx):
    """Return (natoms, info, atom_lines) for the requested frame index."""
    pos = 0
    for _ in range(frame_idx + 1):
        if pos >= len(lines):
            raise IndexError(f"Frame {frame_idx} not found.")
        natoms = int(lines[pos].strip())
        info = _parse_comment(lines[pos + 1])
        atom_lines = lines[pos + 2: pos + 2 + natoms]
        pos += 2 + natoms
    return natoms, info, atom_lines


def parse_atoms(atom_lines, props):
    """Return dict: id -> (x, y, z)  (and also list for ordered access)."""
    col, col_map = 0, {}
    for name, dtype, ncols in props:
        col_map[name] = (col, dtype, ncols)
        col += ncols

    atoms = {}   # id -> (x, y, z)
    order = []   # ids in file order

    for line in atom_lines:
        v = line.split()
        # position
        c, _, n = col_map.get('pos', (None, None, None))
        if c is None:
            raise ValueError("No 'pos' column in properties.")
        xyz = (float(v[c]), float(v[c+1]), float(v[c+2]))
        # atom id
        if 'id' in col_map:
            c_id, _, _ = col_map['id']
            aid = int(v[c_id])
        else:
            aid = len(order) + 1
        atoms[aid] = xyz
        order.append(aid)

    return atoms, order


# ── Geometry ─────────────────────────────────────────────────────────────────

def lattice_to_cell(lattice):
    """Return (a_vec, b_vec, c_vec) as tuples from the 9-value Lattice."""
    return tuple(lattice[0:3]), tuple(lattice[3:6]), tuple(lattice[6:9])


def min_image(dx, dy, dz, a, b, c):
    """Apply minimum image convention for a general triclinic cell.

    Uses fractional coordinates: reduce displacement to [-0.5, 0.5) in
    fractional space, then convert back.
    """
    # Build cell matrix (rows = lattice vectors)
    # Solve s = H^-1 d  (fractional components)
    # For orthorhombic this is trivial; use the general formula otherwise.
    ax, ay, az = a
    bx, by, bz = b
    cx, cy, cz = c

    # Volume (scalar triple product)
    vol = (ax * (by*cz - bz*cy)
         - ay * (bx*cz - bz*cx)
         + az * (bx*cy - by*cx))

    if abs(vol) < 1e-20:
        return dx, dy, dz   # degenerate cell, skip

    # Inverse of the cell matrix (transposed cofactors / vol)
    inv = [
        [(by*cz - bz*cy)/vol, (az*cy - ay*cz)/vol, (ay*bz - az*by)/vol],
        [(bz*cx - bx*cz)/vol, (ax*cz - az*cx)/vol, (az*bx - ax*bz)/vol],
        [(bx*cy - by*cx)/vol, (ay*cx - ax*cy)/vol, (ax*by - ay*bx)/vol],
    ]

    # Fractional displacement
    sa = inv[0][0]*dx + inv[0][1]*dy + inv[0][2]*dz
    sb = inv[1][0]*dx + inv[1][1]*dy + inv[1][2]*dz
    sc = inv[2][0]*dx + inv[2][1]*dy + inv[2][2]*dz

    # Wrap to [-0.5, 0.5)
    sa -= round(sa)
    sb -= round(sb)
    sc -= round(sc)

    # Back to Cartesian
    rx = sa*ax + sb*bx + sc*cx
    ry = sa*ay + sb*by + sc*cy
    rz = sa*az + sb*bz + sc*cz
    return rx, ry, rz


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compute RMSD between two extended XYZ configurations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('file_a', help="First XYZ file")
    parser.add_argument('file_b', help="Second XYZ file")
    parser.add_argument('--frame-a', type=int, default=0,
                        help="Frame index in file A (default: 0)")
    parser.add_argument('--frame-b', type=int, default=0,
                        help="Frame index in file B (default: 0)")
    parser.add_argument('--no-pbc', action='store_true',
                        help="Disable periodic boundary conditions")
    parser.add_argument('--verbose', '-v', action='store_true',
                        help="Print per-atom displacements")
    parser.add_argument('--top', type=int, default=0, metavar='N',
                        help="Print the N atoms with largest displacement")
    args = parser.parse_args()

    # Read files
    try:
        lines_a = open(args.file_a).readlines()
        lines_b = open(args.file_b).readlines()
    except FileNotFoundError as e:
        sys.exit(f"Error: {e}")

    try:
        natoms_a, info_a, atom_lines_a = read_frame(lines_a, args.frame_a)
        natoms_b, info_b, atom_lines_b = read_frame(lines_b, args.frame_b)
    except (IndexError, ValueError) as e:
        sys.exit(f"Error reading frame: {e}")

    props_a = info_a.get('properties')
    props_b = info_b.get('properties')
    if not props_a or not props_b:
        sys.exit("Error: 'properties' key missing from one of the XYZ headers.")

    atoms_a, order_a = parse_atoms(atom_lines_a, props_a)
    atoms_b, order_b = parse_atoms(atom_lines_b, props_b)

    # Match by id
    common_ids = sorted(set(atoms_a) & set(atoms_b))
    if not common_ids:
        sys.exit("Error: no matching atom ids between the two files.")
    if len(common_ids) != natoms_a or len(common_ids) != natoms_b:
        print(f"Warning: {natoms_a} atoms in A, {natoms_b} in B, "
              f"{len(common_ids)} matched by id.", file=sys.stderr)

    n = len(common_ids)

    # Cell from file A
    lattice = info_a.get('lattice', [1,0,0, 0,1,0, 0,0,1])
    a_vec, b_vec, c_vec = lattice_to_cell(lattice)
    use_pbc = not args.no_pbc

    # Per-atom displacements
    disp = []   # (|d|, dx, dy, dz, aid)
    for aid in common_ids:
        xa, ya, za = atoms_a[aid]
        xb, yb, zb = atoms_b[aid]
        dx, dy, dz = xb - xa, yb - ya, zb - za
        if use_pbc:
            dx, dy, dz = min_image(dx, dy, dz, a_vec, b_vec, c_vec)
        d = math.sqrt(dx*dx + dy*dy + dz*dz)
        disp.append((d, dx, dy, dz, aid))

    # Statistics
    rmsd = math.sqrt(sum(d**2 for d, *_ in disp) / n)
    rms_x = math.sqrt(sum(dx**2 for _, dx, _, _, _ in disp) / n)
    rms_y = math.sqrt(sum(dy**2 for _, _, dy, _, _ in disp) / n)
    rms_z = math.sqrt(sum(dz**2 for _, _, _, dz, _ in disp) / n)
    d_max, dx_max, dy_max, dz_max, id_max = max(disp, key=lambda x: x[0])
    d_mean = sum(d for d, *_ in disp) / n

    # ── Output ──
    import os
    na = f"{os.path.basename(args.file_a)}[{args.frame_a}]"
    nb = f"{os.path.basename(args.file_b)}[{args.frame_b}]"
    print(f"\nDistance: {na}  →  {nb}")
    print(f"  Atoms matched  : {n}")
    if use_pbc:
        lx = math.sqrt(sum(x**2 for x in a_vec))
        ly = math.sqrt(sum(x**2 for x in b_vec))
        lz = math.sqrt(sum(x**2 for x in c_vec))
        print(f"  Cell (Å)       : {lx:.6f}  {ly:.6f}  {lz:.6f}  [PBC on]")
    else:
        print(f"  PBC            : off")
    print()
    print(f"  RMSD           : {rmsd:.6f} Å")
    print(f"  Mean |disp|    : {d_mean:.6f} Å")
    print(f"  Max  |disp|    : {d_max:.6f} Å  (atom id {id_max},"
          f"  Δ=({dx_max:.4f}, {dy_max:.4f}, {dz_max:.4f}))")
    print(f"  RMS  Δx        : {rms_x:.6f} Å")
    print(f"  RMS  Δy        : {rms_y:.6f} Å")
    print(f"  RMS  Δz        : {rms_z:.6f} Å")

    if 'energy' in info_a and 'energy' in info_b:
        dE = info_b['energy'] - info_a['energy']
        print(f"\n  E(A)           : {info_a['energy']:.10f} eV")
        print(f"  E(B)           : {info_b['energy']:.10f} eV")
        print(f"  ΔE (B - A)     : {dE:+.10f} eV")

    # Optional: top-N displaced atoms
    if args.top > 0:
        sorted_disp = sorted(disp, key=lambda x: -x[0])
        print(f"\n  Top {args.top} displacements:")
        print(f"  {'id':>8}  {'|d| (Å)':>12}  {'Δx':>10}  {'Δy':>10}  {'Δz':>10}")
        for d, dx, dy, dz, aid in sorted_disp[:args.top]:
            print(f"  {aid:>8}  {d:>12.6f}  {dx:>10.6f}  {dy:>10.6f}  {dz:>10.6f}")

    # Optional: all atoms
    if args.verbose:
        print(f"\n  {'id':>8}  {'|d| (Å)':>12}  {'Δx':>10}  {'Δy':>10}  {'Δz':>10}")
        for d, dx, dy, dz, aid in disp:
            print(f"  {aid:>8}  {d:>12.6f}  {dx:>10.6f}  {dy:>10.6f}  {dz:>10.6f}")

    print()


if __name__ == '__main__':
    main()
