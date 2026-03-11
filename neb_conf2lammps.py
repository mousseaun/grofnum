#!/usr/bin/env python3
"""
Reformat a NEB final-configuration file into LAMMPS data format.

The NEB final file (neb_fin_*.conf) contains atom positions without box
information or atom IDs.  The companion NEB initial file (neb_ini_*.conf)
is a proper LAMMPS data file from which the box bounds, atom IDs, and atom
types are taken.  Atoms in both files are assumed to be in the same order.

Output: a complete LAMMPS data file combining the metadata from the initial
file with the positions from the final file.

Usage:
  ./neb_conf2lammps.py  neb_ini.conf  neb_fin.conf
  ./neb_conf2lammps.py  neb_ini.conf  neb_fin.conf  -o output.data
"""

import argparse
import re
import sys
import os


# ── Parsers ───────────────────────────────────────────────────────────────────

def parse_ini(path):
    """Parse a NEB initial (LAMMPS-format) conf file.

    Returns:
        header_lines : lines before the Atoms section (list of str, no newline)
        atoms        : list of (atom_id, atom_type) in file order
    """
    header_lines = []
    atoms = []
    section = None

    with open(path) as f:
        lines = f.readlines()

    for line in lines:
        s = line.rstrip('\n')
        stripped = s.strip()

        if stripped.lower() == 'atoms':
            section = 'atoms'
            # Don't add "Atoms" to header_lines; write_lammps emits it
            continue

        if section != 'atoms':
            header_lines.append(s)
        else:
            if not stripped or stripped.startswith('#'):
                header_lines.append(s)   # blank/comment after "Atoms" header
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    aid  = int(parts[0])
                    atyp = int(parts[1])
                    atoms.append((aid, atyp))
                    continue
                except ValueError:
                    pass
            # Non-atom line inside Atoms section (shouldn't happen, keep it)
            header_lines.append(s)

    # Trim trailing blank lines from header (they'll be re-added before Atoms)
    while header_lines and not header_lines[-1].strip():
        header_lines.pop()

    return header_lines, atoms


def parse_fin(path):
    """Parse a NEB final conf file: first line = natoms, rest = type x y z.

    Returns list of (atom_type, x, y, z) in file order.
    """
    atoms = []
    with open(path) as f:
        lines = f.readlines()

    if not lines:
        raise ValueError(f"{path!r} is empty.")

    try:
        natoms = int(lines[0].strip())
    except ValueError:
        raise ValueError(
            f"First line of {path!r} should be the atom count, "
            f"got: {lines[0].strip()!r}"
        )

    for i, line in enumerate(lines[1:], start=2):
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        parts = stripped.split()
        if len(parts) < 4:
            raise ValueError(
                f"{path!r} line {i}: expected 'type x y z', got {stripped!r}"
            )
        try:
            atyp = int(parts[0])
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append((atyp, x, y, z))
        except ValueError as e:
            raise ValueError(f"{path!r} line {i}: {e}")

    if len(atoms) != natoms:
        raise ValueError(
            f"{path!r}: header says {natoms} atoms but found {len(atoms)}."
        )
    return atoms


# ── Writer ────────────────────────────────────────────────────────────────────

def write_lammps(out_path, header_lines, ini_atoms, fin_atoms, src_ini, src_fin):
    """Write the combined LAMMPS data file."""
    if len(ini_atoms) != len(fin_atoms):
        raise ValueError(
            f"Atom count mismatch: ini has {len(ini_atoms)}, "
            f"fin has {len(fin_atoms)}."
        )

    with open(out_path, 'w') as f:
        # Header: replace first comment line to document the source
        for i, line in enumerate(header_lines):
            if i == 0:
                # Keep original comment but note the final positions
                orig = line.strip().lstrip('#').strip()
                f.write(f"# {orig}  [final positions from {os.path.basename(src_fin)}]\n")
            else:
                f.write(line + '\n')

        # Atoms section
        f.write('\n Atoms\n\n')
        for (aid, _atyp_ini), (atyp_fin, x, y, z) in zip(ini_atoms, fin_atoms):
            f.write(f"{aid:6d} {atyp_fin:2d}  {x:18.10f}  {y:18.10f}  {z:18.10f}\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Convert NEB final conf to LAMMPS data format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('ini_file',
                        help="NEB initial conf (LAMMPS format, provides box & atom IDs)")
    parser.add_argument('fin_file',
                        help="NEB final conf (type x y z per atom)")
    parser.add_argument('-o', '--output', default=None,
                        help="Output file (default: fin_file with .data extension)")
    args = parser.parse_args()

    # Default output name
    if args.output:
        out_path = args.output
    else:
        base = os.path.splitext(args.fin_file)[0]
        out_path = base + '.data'

    # Parse
    try:
        header_lines, ini_atoms = parse_ini(args.ini_file)
    except (FileNotFoundError, ValueError) as e:
        sys.exit(f"Error reading ini file: {e}")

    try:
        fin_atoms = parse_fin(args.fin_file)
    except (FileNotFoundError, ValueError) as e:
        sys.exit(f"Error reading fin file: {e}")

    if len(ini_atoms) != len(fin_atoms):
        sys.exit(
            f"Error: ini has {len(ini_atoms)} atoms, "
            f"fin has {len(fin_atoms)} atoms – files must have the same atom count."
        )

    # Write
    try:
        write_lammps(out_path, header_lines, ini_atoms, fin_atoms,
                     args.ini_file, args.fin_file)
    except (OSError, ValueError) as e:
        sys.exit(f"Error writing output: {e}")

    # Summary
    ntypes = len({t for _, t in ini_atoms})
    print(f"Wrote: {out_path}")
    print(f"  Atoms     : {len(ini_atoms)}")
    print(f"  Atom types: {ntypes}")
    print(f"  Box from  : {os.path.basename(args.ini_file)}")
    print(f"  Positions : {os.path.basename(args.fin_file)}")


if __name__ == '__main__':
    main()
