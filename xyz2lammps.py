#!/usr/bin/env python3
"""
Convert extended XYZ file to LAMMPS data file format.

Extended XYZ format (line 2 example):
  Lattice="ax ay az bx by bz cx cy cz" properties=species:S:1:pos:R:3:forces:R:3:id:I:1 energy=-1234.5

Usage:
  ./xyz2lammps.py input.xyz                              # -> input.data
  ./xyz2lammps.py input.xyz -o output.data
  ./xyz2lammps.py traj.xyz --frame 5 -o frame5.data
  ./xyz2lammps.py input.xyz --symbols Si                 # one type -> Si
  ./xyz2lammps.py input.xyz --symbols Al,O               # two types -> Al=1, O=2
  ./xyz2lammps.py traj.xyz --list-frames
"""

import argparse
import re
import sys
import math
from collections import OrderedDict

# -----------------------------------------------------------------------------
# Atomic masses (g/mol)
# -----------------------------------------------------------------------------
ATOMIC_MASSES = {
    'H': 1.008,    'He': 4.003,   'Li': 6.941,   'Be': 9.012,
    'B': 10.811,   'C': 12.011,   'N': 14.007,   'O': 15.999,
    'F': 18.998,   'Ne': 20.180,  'Na': 22.990,  'Mg': 24.305,
    'Al': 26.982,  'Si': 28.086,  'P': 30.974,   'S': 32.065,
    'Cl': 35.453,  'Ar': 39.948,  'K': 39.098,   'Ca': 40.078,
    'Sc': 44.956,  'Ti': 47.867,  'V': 50.942,   'Cr': 51.996,
    'Mn': 54.938,  'Fe': 55.845,  'Co': 58.933,  'Ni': 58.693,
    'Cu': 63.546,  'Zn': 65.380,  'Ga': 69.723,  'Ge': 72.630,
    'As': 74.922,  'Se': 78.971,  'Br': 79.904,  'Kr': 83.798,
    'Rb': 85.468,  'Sr': 87.620,  'Y': 88.906,   'Zr': 91.224,
    'Nb': 92.906,  'Mo': 95.960,  'Tc': 98.000,  'Ru': 101.070,
    'Rh': 102.906, 'Pd': 106.420, 'Ag': 107.868, 'Cd': 112.411,
    'In': 114.818, 'Sn': 118.710, 'Sb': 121.760, 'Te': 127.600,
    'I': 126.904,  'Xe': 131.293, 'Cs': 132.905, 'Ba': 137.327,
    'La': 138.905, 'Ce': 140.116, 'Hf': 178.490, 'Ta': 180.948,
    'W': 183.840,  'Re': 186.207, 'Os': 190.230, 'Ir': 192.217,
    'Pt': 195.084, 'Au': 196.967, 'Hg': 200.592, 'Tl': 204.383,
    'Pb': 207.200, 'Bi': 208.980,
}

# -----------------------------------------------------------------------------
# Parsing
# -----------------------------------------------------------------------------

def parse_properties(prop_str):
    """Parse 'name:type:ncols:...' into list of (name, type, ncols)."""
    parts = prop_str.strip().split(':')
    if len(parts) % 3 != 0:
        raise ValueError(f"Malformed properties string: {prop_str!r}")
    result = []
    for i in range(0, len(parts), 3):
        name = parts[i]
        dtype = parts[i + 1].upper()
        ncols = int(parts[i + 2])
        result.append((name.lower(), dtype, ncols))
    return result


def parse_comment_line(line):
    """Parse extended XYZ comment (header) line.

    Returns dict with keys: 'lattice' (list of 9 floats), 'properties'
    (list of tuples), 'energy' (float, optional).
    """
    info = {}

    # Lattice
    m = re.search(r'[Ll]attice\s*=\s*"([^"]*)"', line)
    if m:
        info['lattice'] = list(map(float, m.group(1).split()))

    # Properties
    m = re.search(r'[Pp]roperties\s*=\s*(\S+)', line)
    if m:
        info['properties'] = parse_properties(m.group(1))

    # Energy  (energy=value  or  energy:value)
    m = re.search(r'[Ee]nergy\s*[=:]\s*(-?[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)', line)
    if m:
        info['energy'] = float(m.group(1))

    return info


def read_frame(lines, start):
    """Read one XYZ frame starting at line index `start`.

    Returns (natoms, comment_info, atom_lines, next_start).
    """
    natoms = int(lines[start].strip())
    info = parse_comment_line(lines[start + 1])
    atom_lines = lines[start + 2: start + 2 + natoms]
    return natoms, info, atom_lines, start + 2 + natoms


def parse_atoms(atom_lines, props):
    """Parse atom data lines according to properties descriptor.

    Returns dict with keys 'species', 'pos', and optionally 'forces', 'id'.
    """
    # Build column index map: name -> starting column
    col_start = {}
    col = 0
    for name, dtype, ncols in props:
        col_start[name] = (col, dtype, ncols)
        col += ncols

    data = {'species': [], 'pos': []}

    for line in atom_lines:
        vals = line.split()

        # Species
        if 'species' in col_start:
            c, dtype, n = col_start['species']
            data['species'].append(vals[c])

        # Position
        if 'pos' in col_start:
            c, dtype, n = col_start['pos']
            data['pos'].append([float(vals[c + i]) for i in range(3)])

        # Forces (optional)
        for fname in ('forces', 'force'):
            if fname in col_start:
                if 'forces' not in data:
                    data['forces'] = []
                c, dtype, n = col_start[fname]
                data['forces'].append([float(vals[c + i]) for i in range(3)])
                break

        # Atom id (optional)
        if 'id' in col_start:
            if 'id' not in data:
                data['id'] = []
            c, dtype, n = col_start['id']
            data['id'].append(int(vals[c]))

    return data


# -----------------------------------------------------------------------------
# Box conversion: extended XYZ lattice -> LAMMPS triclinic parameters
# -----------------------------------------------------------------------------

def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def norm(v):
    return math.sqrt(dot(v, v))


def lattice_to_lammps(lattice):
    """Convert the 9-value Lattice (row vectors a, b, c) to LAMMPS
    triclinic box parameters (lx, ly, lz, xy, xz, yz).

    LAMMPS requires the upper-triangular form:
        a = (lx, 0,  0 )
        b = (xy, ly, 0 )
        c = (xz, yz, lz)
    """
    a = lattice[0:3]
    b = lattice[3:6]
    c = lattice[6:9]

    lx = norm(a)
    if lx < 1e-12:
        raise ValueError("Zero-length lattice vector a")

    a_hat = [x / lx for x in a]

    xy = dot(b, a_hat)
    b_perp = [b[i] - xy * a_hat[i] for i in range(3)]
    ly = norm(b_perp)

    xz = dot(c, a_hat)
    if ly > 1e-12:
        b_perp_hat = [x / ly for x in b_perp]
        yz = dot(c, b_perp_hat)
    else:
        yz = 0.0

    lz_sq = dot(c, c) - xz**2 - yz**2
    lz = math.sqrt(max(0.0, lz_sq))

    return lx, ly, lz, xy, xz, yz


# -----------------------------------------------------------------------------
# Species mapping
# -----------------------------------------------------------------------------

def build_type_map(species_list, props, user_symbols=None):
    """Map species values to LAMMPS integer types (starting from 1).

    If species are stored as integers in the XYZ file, the unique values are
    sorted and renumbered 1..N.  If they are element symbols (strings), they
    are sorted alphabetically and renumbered.

    When `user_symbols` is provided (list of element symbols like ['Si'] or
    ['Al', 'O']), the type names used for the Masses section are taken from
    there instead of the raw species values.  This is useful when the XYZ file
    stores integer type indices (1, 2, ...) instead of element symbols.

    Returns:
        type_map  : dict  species_str -> lammps_type (int)
        type_names: list  element symbol (or raw species) per type, 0-indexed
    """
    # Determine the dtype of species from properties descriptor
    species_dtype = 'S'
    for name, dtype, ncols in props:
        if name == 'species':
            species_dtype = dtype
            break

    unique = list(OrderedDict.fromkeys(species_list))

    if species_dtype == 'I':
        unique_sorted = sorted(unique, key=lambda x: int(x))
    else:
        unique_sorted = sorted(unique)

    type_map = {sp: i + 1 for i, sp in enumerate(unique_sorted)}

    # Override type names with user-supplied element symbols
    if user_symbols:
        if len(user_symbols) != len(unique_sorted):
            raise ValueError(
                f"--symbols has {len(user_symbols)} entries but {len(unique_sorted)} "
                f"distinct species found in file: {unique_sorted}"
            )
        type_names = user_symbols
    else:
        type_names = unique_sorted

    return type_map, type_names


# -----------------------------------------------------------------------------
# Writer
# -----------------------------------------------------------------------------

def write_lammps_data(out_path, natoms, info, atom_data, type_map,
                      type_names, lx, ly, lz, xy, xz, yz, src_name):
    """Write LAMMPS data file in atom_style atomic."""
    ntypes = len(type_map)
    is_triclinic = (abs(xy) > 1e-10 or abs(xz) > 1e-10 or abs(yz) > 1e-10)

    atom_ids = atom_data.get('id') or list(range(1, natoms + 1))

    with open(out_path, 'w') as f:
        # Header comment
        energy_str = (f"  energy = {info['energy']:.10f}"
                      if 'energy' in info else "")
        f.write(f"LAMMPS data file – converted from {src_name}{energy_str}\n\n")

        f.write(f"{natoms} atoms\n")
        f.write(f"{ntypes} atom types\n\n")

        # Box
        f.write(f"0.0 {lx:.10f}  xlo xhi\n")
        f.write(f"0.0 {ly:.10f}  ylo yhi\n")
        f.write(f"0.0 {lz:.10f}  zlo zhi\n")
        if is_triclinic:
            f.write(f"{xy:.10f} {xz:.10f} {yz:.10f}  xy xz yz\n")

        # Masses
        f.write("\nMasses\n\n")
        for t, sp in enumerate(type_names, start=1):
            mass = ATOMIC_MASSES.get(sp, 0.0)
            comment = f"  # {sp}"
            if mass == 0.0:
                comment += "  (unknown – please set manually)"
            f.write(f"{t} {mass:.4f}{comment}\n")

        # Atoms
        f.write("\nAtoms  # atomic\n\n")
        for i in range(natoms):
            aid = atom_ids[i]
            sp = atom_data['species'][i]
            t = type_map[sp]
            x, y_coord, z = atom_data['pos'][i]
            f.write(f"{aid:8d} {t} {x:18.10f} {y_coord:18.10f} {z:18.10f}\n")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert extended XYZ file to LAMMPS data file (atom_style atomic).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('input', help="Input extended XYZ file")
    parser.add_argument('-o', '--output', default=None,
                        help="Output LAMMPS data file (default: <input>.data)")
    parser.add_argument('--frame', type=int, default=0,
                        help="Frame index to convert (0-based, default: 0)")
    parser.add_argument('--list-frames', action='store_true',
                        help="List all frames found in the file and exit")
    parser.add_argument('--symbols', default=None,
                        help=("Comma-separated element symbols in type order "
                              "(e.g. 'Si' or 'Al,O'). Use when the XYZ file "
                              "stores integer species indices instead of element names."))
    args = parser.parse_args()

    # Default output name
    if args.output is None:
        base = args.input
        if base.endswith('.xyz'):
            base = base[:-4]
        args.output = base + '.data'

    # Read file
    try:
        with open(args.input, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        sys.exit(f"Error: file not found: {args.input}")

    # Scan frames
    if args.list_frames:
        idx = 0
        frame_idx = 0
        while idx < len(lines):
            try:
                natoms = int(lines[idx].strip())
            except (ValueError, IndexError):
                break
            info = parse_comment_line(lines[idx + 1])
            energy_str = (f"  energy = {info['energy']:.6f}"
                          if 'energy' in info else "")
            print(f"Frame {frame_idx:4d}: {natoms:6d} atoms{energy_str}")
            idx += 2 + natoms
            frame_idx += 1
        return

    # Read the requested frame
    idx = 0
    for current in range(args.frame + 1):
        if idx >= len(lines):
            sys.exit(f"Error: frame {args.frame} not found (file has {current} frames).")
        try:
            natoms, info, atom_lines, idx = read_frame(lines, idx)
        except (ValueError, IndexError) as e:
            sys.exit(f"Error reading frame {current}: {e}")

    # Check properties
    props = info.get('properties')
    if props is None:
        sys.exit("Error: no 'properties' key found in XYZ header.")

    # Parse atoms
    atom_data = parse_atoms(atom_lines, props)
    if not atom_data['species']:
        sys.exit("Error: no species column found in atom data.")
    if not atom_data['pos']:
        sys.exit("Error: no pos column found in atom data.")

    # Lattice
    lattice = info.get('lattice')
    if lattice is None or len(lattice) != 9:
        sys.exit("Error: Lattice not found or malformed in XYZ header.")
    lx, ly, lz, xy, xz, yz = lattice_to_lammps(lattice)

    # Type mapping
    user_symbols = [s.strip() for s in args.symbols.split(',')] if args.symbols else None
    try:
        type_map, type_names = build_type_map(atom_data['species'], props, user_symbols)
    except ValueError as e:
        sys.exit(f"Error: {e}")

    # Write
    import os
    write_lammps_data(
        args.output, natoms, info, atom_data, type_map, type_names,
        lx, ly, lz, xy, xz, yz,
        src_name=os.path.basename(args.input),
    )

    # Summary
    ntypes = len(type_map)
    is_triclinic = (abs(xy) > 1e-10 or abs(xz) > 1e-10 or abs(yz) > 1e-10)
    box_type = "triclinic" if is_triclinic else "orthorhombic"
    type_summary = "  ".join(f"{sp}→{t}" for sp, t in
                             sorted(type_map.items(), key=lambda x: x[1]))
    print(f"Wrote: {args.output}")
    print(f"  atoms   : {natoms}")
    print(f"  types   : {ntypes}  [{type_summary}]")
    print(f"  box     : {box_type}  lx={lx:.6f}  ly={ly:.6f}  lz={lz:.6f}")
    if is_triclinic:
        print(f"  tilt    : xy={xy:.6f}  xz={xz:.6f}  yz={yz:.6f}")
    if 'energy' in info:
        print(f"  energy  : {info['energy']:.10f}")
    if 'forces' in atom_data:
        print("  (forces from XYZ not written – not part of LAMMPS data format)")


if __name__ == '__main__':
    main()
