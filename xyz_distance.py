#!/usr/bin/env python3
"""
Compute the structural distance between two configurations.

Supported formats: extended XYZ and LAMMPS data (atom_style atomic).
Format is auto-detected from the file extension (.xyz → XYZ; .data /
.lammps / .lmp → LAMMPS) and can be overridden with --format-a / --format-b.

The distance is the RMSD (root-mean-square displacement) between matching
atoms.  Atoms are matched by their 'id' field when present, otherwise by
line order.  Periodic boundary conditions (minimum image convention) are
applied using the lattice from file A.

Outputs:
  - RMSD (Å)
  - Max per-atom displacement and which atom
  - RMS per component (dx, dy, dz)
  - Energy difference (if both files carry an energy)

For XYZ trajectory files, use --frame-a / --frame-b to select frames.

Usage:
  ./xyz_distance.py A.xyz B.xyz
  ./xyz_distance.py A.data B.data
  ./xyz_distance.py A.xyz B.data
  ./xyz_distance.py traj.xyz B.xyz --frame-a 5
  ./xyz_distance.py A.xyz B.xyz --no-pbc
  ./xyz_distance.py A.xyz B.xyz --verbose       # print per-atom displacements
  ./xyz_distance.py A.cfg B.cfg --format-a lammps --format-b lammps
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


# ── LAMMPS parser ────────────────────────────────────────────────────────────

def read_lammps(path):
    """Parse a LAMMPS data file (atom_style atomic).

    Returns (natoms, info, atoms, order) where:
      info  : dict with 'lattice' (9 floats, row-major) and optionally 'energy'
      atoms : dict  id -> (x, y, z)
      order : list of ids in file order
    """
    with open(path) as f:
        lines = f.readlines()

    info = {}
    natoms = None
    xlo = xhi = ylo = yhi = zlo = zhi = 0.0
    xy = xz = yz = 0.0
    atoms = {}
    order = []
    section = None

    # Energy may appear in the header comment (produced by xyz2lammps.py)
    if lines:
        m = re.search(
            r'energy\s*=\s*(-?[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)',
            lines[0]
        )
        if m:
            info['energy'] = float(m.group(1))

    for line in lines:
        stripped = line.strip()
        # Strip inline comments for keyword lines, but keep raw for atom data
        no_comment = stripped.split('#')[0].strip()

        if not stripped or stripped.startswith('#'):
            continue

        low = stripped.lower()

        # Section headers (must be checked before numeric parsing)
        if low.startswith('atoms'):
            section = 'atoms'
            continue
        if low == 'masses':
            section = 'masses'
            continue
        if low in ('velocities', 'bonds', 'angles', 'dihedrals', 'impropers',
                   'pair coeffs', 'bond coeffs', 'angle coeffs'):
            section = low
            continue

        # Skip non-Atoms sections
        if section != 'atoms':
            # Atom count
            m = re.match(r'^(\d+)\s+atoms\s*$', no_comment)
            if m:
                natoms = int(m.group(1))
                continue

            # Box bounds
            m = re.match(
                r'^(-?[\d.eE+\-]+)\s+(-?[\d.eE+\-]+)\s+xlo\s+xhi', no_comment)
            if m:
                xlo, xhi = float(m.group(1)), float(m.group(2))
                continue
            m = re.match(
                r'^(-?[\d.eE+\-]+)\s+(-?[\d.eE+\-]+)\s+ylo\s+yhi', no_comment)
            if m:
                ylo, yhi = float(m.group(1)), float(m.group(2))
                continue
            m = re.match(
                r'^(-?[\d.eE+\-]+)\s+(-?[\d.eE+\-]+)\s+zlo\s+zhi', no_comment)
            if m:
                zlo, zhi = float(m.group(1)), float(m.group(2))
                continue
            m = re.match(
                r'^(-?[\d.eE+\-]+)\s+(-?[\d.eE+\-]+)\s+(-?[\d.eE+\-]+)'
                r'\s+xy\s+xz\s+yz', no_comment)
            if m:
                xy, xz, yz = float(m.group(1)), float(m.group(2)), float(m.group(3))
                continue
        else:
            # Atom lines: atom_id  atom_type  x  y  z  [ix iy iz]
            parts = stripped.split()
            if len(parts) >= 5:
                try:
                    aid = int(parts[0])
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atoms[aid] = (x, y, z)
                    order.append(aid)
                except ValueError:
                    pass  # skip malformed lines

    # Build 9-value lattice from LAMMPS triclinic box:
    #   a = (lx,  0,  0)
    #   b = (xy, ly,  0)
    #   c = (xz, yz, lz)
    lx, ly, lz = xhi - xlo, yhi - ylo, zhi - zlo
    info['lattice'] = [lx, 0.0, 0.0,
                       xy,  ly, 0.0,
                       xz,  yz,  lz]

    if natoms is None:
        natoms = len(atoms)

    return natoms, info, atoms, order


# ── Format detection & unified loader ────────────────────────────────────────

_LAMMPS_EXTS = {'.data', '.lammps', '.lmp'}
_XYZ_EXTS    = {'.xyz'}


def _detect_format(path, override):
    if override:
        return override
    import os
    ext = os.path.splitext(path)[1].lower()
    if ext in _LAMMPS_EXTS:
        return 'lammps'
    return 'xyz'


def load_frame(path, fmt, frame_idx):
    """Load one frame from *path* and return (natoms, info, atoms, order).

    atoms : dict  id -> (x, y, z)
    order : list of ids in file order
    """
    if fmt == 'lammps':
        if frame_idx != 0:
            print("Warning: LAMMPS data files are single-frame; "
                  "--frame option ignored.", file=sys.stderr)
        return read_lammps(path)

    # XYZ
    lines = open(path).readlines()
    natoms, info, atom_lines = read_frame(lines, frame_idx)
    props = info.get('properties')
    if not props:
        sys.exit(f"Error: 'properties' key missing from XYZ header in {path}.")
    atoms, order = parse_atoms(atom_lines, props)
    return natoms, info, atoms, order


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
        description="Compute RMSD between two atomic configurations (XYZ or LAMMPS).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('file_a', help="First file (XYZ or LAMMPS data)")
    parser.add_argument('file_b', help="Second file (XYZ or LAMMPS data)")
    parser.add_argument('--frame-a', type=int, default=0,
                        help="Frame index in file A, XYZ only (default: 0)")
    parser.add_argument('--frame-b', type=int, default=0,
                        help="Frame index in file B, XYZ only (default: 0)")
    parser.add_argument('--format-a', choices=['xyz', 'lammps'], default=None,
                        help="Force format for file A (default: auto from extension)")
    parser.add_argument('--format-b', choices=['xyz', 'lammps'], default=None,
                        help="Force format for file B (default: auto from extension)")
    parser.add_argument('--no-pbc', action='store_true',
                        help="Disable periodic boundary conditions")
    parser.add_argument('--verbose', '-v', action='store_true',
                        help="Print per-atom displacements")
    parser.add_argument('--top', type=int, default=0, metavar='N',
                        help="Print the N atoms with largest displacement")
    args = parser.parse_args()

    fmt_a = _detect_format(args.file_a, args.format_a)
    fmt_b = _detect_format(args.file_b, args.format_b)

    # Read files
    try:
        natoms_a, info_a, atoms_a, order_a = load_frame(
            args.file_a, fmt_a, args.frame_a)
        natoms_b, info_b, atoms_b, order_b = load_frame(
            args.file_b, fmt_b, args.frame_b)
    except FileNotFoundError as e:
        sys.exit(f"Error: {e}")
    except (IndexError, ValueError) as e:
        sys.exit(f"Error reading frame: {e}")

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
    total = math.sqrt(sum(d**2 for d, *_ in disp) )
    rmsd = math.sqrt(sum(d**2 for d, *_ in disp) / n)
    rms_x = math.sqrt(sum(dx**2 for _, dx, _, _, _ in disp) / n)
    rms_y = math.sqrt(sum(dy**2 for _, _, dy, _, _ in disp) / n)
    rms_z = math.sqrt(sum(dz**2 for _, _, _, dz, _ in disp) / n)
    d_max, dx_max, dy_max, dz_max, id_max = max(disp, key=lambda x: x[0])
    d_mean = sum(d for d, *_ in disp) / n

    # ── Output ──
    import os
    def _label(path, fmt, frame):
        base = os.path.basename(path)
        if fmt == 'lammps':
            return f"{base} [lammps]"
        return f"{base}[{frame}]"

    na = _label(args.file_a, fmt_a, args.frame_a)
    nb = _label(args.file_b, fmt_b, args.frame_b)
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
    print(f"  Distance totale   : {total:.6f} Å")
    print(f"  RMSD              : {rmsd:.6f} Å")
    print(f"  Mean |disp|       : {d_mean:.6f} Å")
    print(f"  Max  |disp|       : {d_max:.6f} Å  (atom id {id_max},"
          f"  Δ=({dx_max:.4f}, {dy_max:.4f}, {dz_max:.4f}))")
    print(f"  RMS  Δx           : {rms_x:.6f} Å")
    print(f"  RMS  Δy           : {rms_y:.6f} Å")
    print(f"  RMS  Δz           : {rms_z:.6f} Å")

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
