#!/usr/bin/env python3
"""
Apply an ARTn event to a LAMMPS data file.

Reads three files:
  1. A LAMMPS data file        – initial configuration (all atoms)
  2. An ARTn event file        – cluster atom IDs (ncluster + list)
  3. An event XYZ file         – 3 frames: Initial / Saddle / Final
                                 (positions of the cluster atoms only)

For each cluster atom i (local index in XYZ ↔ cluster_ids[i-1] in LAMMPS):
    disp_saddle[i] = frame2[i] − frame1[i]
    disp_final[i]  = frame3[i] − frame1[i]

The displacements are added to the atom's current position in the LAMMPS
file.  Periodic boundary conditions (minimum image then wrap) are applied
by default.

Outputs (two LAMMPS data files):
  <base>_saddle.data   – initial structure with saddle displacements applied
  <base>_final.data    – initial structure with final  displacements applied

Usage:
  ./apply_event.py init.data event_dir/event event_dir/event.xyz
  ./apply_event.py init.data event_dir/event event_dir/event.xyz -o result
  ./apply_event.py init.data event_dir/event event_dir/event.xyz --no-pbc -v
"""

import argparse
import math
import re
import sys
import os


# ── Event text-file parser ────────────────────────────────────────────────────

def parse_event_file(path):
    """Return (ncluster, cluster_ids).

    cluster_ids[i] is the LAMMPS atom ID corresponding to local XYZ index i+1.
    Zeros in the padded array are ignored; exactly ncluster IDs are returned.
    """
    ncluster = None
    shape    = None
    cluster_ids = []
    in_cluster  = False
    n_read      = 0

    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if in_cluster:
                if re.match(r'^-?\d+$', s):
                    val = int(s)
                    n_read += 1
                    if val != 0 and len(cluster_ids) < (ncluster or 0):
                        cluster_ids.append(val)
                    if shape and n_read >= shape:
                        break
                elif n_read > 0:
                    break  # first non-numeric line ends the section
                continue

            m = re.match(r'ncluster\s*:\s*(\d+)', s)
            if m:
                ncluster = int(m.group(1))
                continue

            m = re.match(r'shape of cluster atoms\s*:\s*(\d+)', s)
            if m:
                shape = int(m.group(1))
                continue

            if re.match(r'cluster atoms\s*:', s) and 'shape' not in s:
                in_cluster = True

    if ncluster is None:
        raise ValueError(f"'ncluster' not found in {path!r}.")
    if len(cluster_ids) < ncluster:
        raise ValueError(
            f"Expected {ncluster} cluster atom IDs, found {len(cluster_ids)} "
            f"in {path!r}."
        )
    return ncluster, cluster_ids[:ncluster]


# ── Event XYZ parser ──────────────────────────────────────────────────────────

def parse_event_xyz(path):
    """Return list of frames.

    Each frame is a dict: local_id (1-based int) -> (x, y, z).
    The comment line is also parsed for ConfE and DelE.
    """
    frames = []
    with open(path) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            i += 1
            continue
        try:
            natoms = int(s)
        except ValueError:
            i += 1
            continue

        comment = lines[i + 1]
        conf_e = del_e = None
        m = re.search(r'ConfE\s*=\s*"?\s*(-?[\d.eE+\-]+)', comment)
        if m:
            conf_e = float(m.group(1))
        m = re.search(r'DelE\s*=\s*"?\s*(-?[\d.eE+\-]+)', comment)
        if m:
            del_e = float(m.group(1))

        atoms = {}
        for j in range(natoms):
            parts = lines[i + 2 + j].split()
            # format: element  x  y  z  local_id
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            lid = int(parts[4])
            atoms[lid] = (x, y, z)

        frames.append({'atoms': atoms, 'conf_e': conf_e, 'del_e': del_e})
        i += 2 + natoms

    return frames


# ── LAMMPS data-file parser ───────────────────────────────────────────────────

def parse_lammps(path):
    """Return (box, atom_positions, raw_lines).

    box            : dict with xlo, xhi, ylo, yhi, zlo, zhi
    atom_positions : dict  atom_id -> (x, y, z)
    raw_lines      : list of original lines (for faithful re-writing)
    """
    box = {}
    atom_positions = {}
    section = None

    with open(path) as f:
        raw_lines = f.readlines()

    for line in raw_lines:
        s = line.strip()
        no_comment = s.split('#')[0].strip()

        if not s or s.startswith('#'):
            continue

        low = s.lower()
        if low.startswith('atoms') and not re.match(r'^\d', s):
            section = 'atoms'
            continue
        if low in ('masses', 'velocities', 'bonds', 'angles',
                   'dihedrals', 'impropers'):
            section = low
            continue

        if section != 'atoms':
            m = re.match(
                r'^(-?[\d.eE+\-]+)\s+(-?[\d.eE+\-]+)\s+(xlo\s+xhi)', no_comment)
            if m:
                box['xlo'], box['xhi'] = float(m.group(1)), float(m.group(2))
                continue
            m = re.match(
                r'^(-?[\d.eE+\-]+)\s+(-?[\d.eE+\-]+)\s+(ylo\s+yhi)', no_comment)
            if m:
                box['ylo'], box['yhi'] = float(m.group(1)), float(m.group(2))
                continue
            m = re.match(
                r'^(-?[\d.eE+\-]+)\s+(-?[\d.eE+\-]+)\s+(zlo\s+zhi)', no_comment)
            if m:
                box['zlo'], box['zhi'] = float(m.group(1)), float(m.group(2))
                continue
        else:
            parts = s.split()
            if len(parts) >= 5:
                try:
                    aid = int(parts[0])
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atom_positions[aid] = (x, y, z)
                except ValueError:
                    pass

    return box, atom_positions, raw_lines


# ── PBC wrapping ──────────────────────────────────────────────────────────────

def wrap(x, lo, hi):
    """Wrap x into [lo, hi) using periodic boundary conditions."""
    L = hi - lo
    return lo + (x - lo) % L


# ── LAMMPS output writer ──────────────────────────────────────────────────────

def write_lammps(out_path, raw_lines, updated_positions, header_comment=None):
    """Write a LAMMPS data file, replacing atom positions from updated_positions.

    updated_positions : dict  atom_id -> (x, y, z)
    Only atoms present in updated_positions are modified; others are unchanged.
    """
    section = None
    with open(out_path, 'w') as f:
        for idx, line in enumerate(raw_lines):
            s = line.strip()
            low = s.lower()

            # Optional header comment replacement (first line only)
            if idx == 0 and header_comment is not None:
                f.write(header_comment + '\n')
                continue

            if low.startswith('atoms') and not re.match(r'^\d', s):
                section = 'atoms'
                f.write(line)
                continue
            if low in ('masses', 'velocities', 'bonds', 'angles',
                       'dihedrals', 'impropers'):
                section = low
                f.write(line)
                continue

            if section == 'atoms' and s:
                parts = s.split()
                if len(parts) >= 5:
                    try:
                        aid = int(parts[0])
                        if aid in updated_positions:
                            nx, ny, nz = updated_positions[aid]
                            # Reconstruct line preserving atom_type (parts[1])
                            f.write(f"{aid:8d} {parts[1]}"
                                    f" {nx:22.16g}"
                                    f" {ny:22.16g}"
                                    f" {nz:22.16g}\n")
                            continue
                    except ValueError:
                        pass
            f.write(line)


# ── Statistics helper ─────────────────────────────────────────────────────────

def displacement_stats(displacements):
    """displacements : list of (lid, dx, dy, dz, |d|)"""
    n = len(displacements)
    rmsd   = math.sqrt(sum(t[4]**2 for t in displacements) / n)
    d_mean = sum(t[4] for t in displacements) / n
    d_max  = max(displacements, key=lambda t: t[4])
    return rmsd, d_mean, d_max


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Apply an ARTn event (Saddle & Final) to a LAMMPS data file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('init_file',  help="Initial LAMMPS data file")
    parser.add_argument('event_file', help="ARTn event text file (cluster atom IDs)")
    parser.add_argument('event_xyz',  help="Event XYZ file (3 frames: Init/Saddle/Final)")
    parser.add_argument('-o', '--output', default=None,
                        help="Base name for output files (default: event file name)")
    parser.add_argument('--saddle-out', default=None,
                        help="Explicit output path for the saddle configuration")
    parser.add_argument('--final-out',  default=None,
                        help="Explicit output path for the final  configuration")
    parser.add_argument('--no-pbc', action='store_true',
                        help="Disable periodic boundary condition wrapping")
    parser.add_argument('--verbose', '-v', action='store_true',
                        help="Print per-atom displacement table")
    parser.add_argument('--top', type=int, default=10, metavar='N',
                        help="Number of largest displacements to show (default: 10)")
    args = parser.parse_args()

    # ── Default output names ──────────────────────────────────────────────────
    base = args.output or os.path.splitext(
        os.path.basename(args.event_file))[0]
    saddle_path = args.saddle_out or f"{base}_saddle.data"
    final_path  = args.final_out  or f"{base}_final.data"

    # ── Parse inputs ──────────────────────────────────────────────────────────
    try:
        ncluster, cluster_ids = parse_event_file(args.event_file)
    except (FileNotFoundError, ValueError) as e:
        sys.exit(f"Error reading event file: {e}")

    try:
        frames = parse_event_xyz(args.event_xyz)
    except FileNotFoundError as e:
        sys.exit(f"Error reading event XYZ: {e}")

    if len(frames) < 3:
        sys.exit(f"Error: expected 3 frames in event XYZ, found {len(frames)}.")

    try:
        box, init_pos, raw_lines = parse_lammps(args.init_file)
    except FileNotFoundError as e:
        sys.exit(f"Error reading LAMMPS file: {e}")

    # ── Build local_id -> LAMMPS atom_id map ─────────────────────────────────
    # local XYZ index is 1-based; cluster_ids is 0-indexed
    local_to_lammps = {i + 1: cluster_ids[i] for i in range(ncluster)}

    # Verify all cluster atoms are present in the LAMMPS file
    missing = [lid for lid in cluster_ids if lid not in init_pos]
    if missing:
        sys.exit(
            f"Error: {len(missing)} cluster atom(s) not found in LAMMPS file: "
            f"{missing[:10]}{'...' if len(missing) > 10 else ''}"
        )

    # ── Compute displacements (frame2−frame1 and frame3−frame1) ──────────────
    f1 = frames[0]['atoms']
    f2 = frames[1]['atoms']
    f3 = frames[2]['atoms']

    disp_saddle = {}   # local_id -> (dx, dy, dz)
    disp_final  = {}

    for lid in range(1, ncluster + 1):
        x1, y1, z1 = f1[lid]
        x2, y2, z2 = f2[lid]
        x3, y3, z3 = f3[lid]
        disp_saddle[lid] = (x2 - x1, y2 - y1, z2 - z1)
        disp_final[lid]  = (x3 - x1, y3 - y1, z3 - z1)

    use_pbc = not args.no_pbc

    # ── Apply displacements to LAMMPS positions ───────────────────────────────
    def apply_displacements(disp):
        updated = {}
        stats   = []
        for lid, lammps_id in local_to_lammps.items():
            x0, y0, z0 = init_pos[lammps_id]
            dx, dy, dz = disp[lid]
            nx, ny, nz = x0 - dx, y0 - dy, z0 - dz
            if use_pbc:
                nx = wrap(nx, box['xlo'], box['xhi'])
                ny = wrap(ny, box['ylo'], box['yhi'])
                nz = wrap(nz, box['zlo'], box['zhi'])
            d = math.sqrt(dx*dx + dy*dy + dz*dz)
            updated[lammps_id] = (nx, ny, nz)
            stats.append((lid, dx, dy, dz, d))
        return updated, stats

    updated_saddle, stats_saddle = apply_displacements(disp_saddle)
    updated_final,  stats_final  = apply_displacements(disp_final)

    # ── Write outputs ─────────────────────────────────────────────────────────
    e1 = frames[0].get('conf_e')
    e2 = frames[1].get('conf_e')
    e3 = frames[2].get('conf_e')

    def _header(label, e):
        estr = f"  energy = {e:.10f}" if e is not None else ""
        src = os.path.basename(args.init_file)
        ev  = os.path.basename(args.event_file)
        return f"LAMMPS data file – {src} + {ev} ({label}){estr}"

    write_lammps(saddle_path, raw_lines, updated_saddle,
                 header_comment=_header("Saddle", e2))
    write_lammps(final_path,  raw_lines, updated_final,
                 header_comment=_header("Final",  e3))

    # ── Summary ───────────────────────────────────────────────────────────────
    rmsd_s, mean_s, max_s = displacement_stats(stats_saddle)
    rmsd_f, mean_f, max_f = displacement_stats(stats_final)

    print(f"\nEvent applied: {os.path.basename(args.event_file)}")
    print(f"  Cluster atoms   : {ncluster}")
    print(f"  PBC wrapping    : {'on' if use_pbc else 'off'}")
    if e1 is not None:
        print(f"  E(Initial)      : {e1:.6f}")

    def _print_config(label, out_path, rmsd, mean, max_t, e):
        print(f"\n  [{label}]  →  {out_path}")
        if e is not None:
            print(f"    ConfE           : {e:.6f}")
        print(f"    RMSD            : {rmsd:.6f} Å")
        print(f"    Mean |disp|     : {mean:.6f} Å")
        print(f"    Max  |disp|     : {max_t[4]:.6f} Å"
              f"  (local atom {max_t[0]}, LAMMPS id {local_to_lammps[max_t[0]]})")

    _print_config("Saddle", saddle_path, rmsd_s, mean_s, max_s, e2)
    _print_config("Final",  final_path,  rmsd_f, mean_f, max_f, e3)

    # Optional: per-atom displacement table
    for label, stats in [("Saddle", stats_saddle), ("Final", stats_final)]:
        top_n = args.top
        sorted_stats = sorted(stats, key=lambda t: -t[4])
        print(f"\n  Top {top_n} displacements [{label}]:")
        print(f"  {'loc':>5}  {'LAMMPS':>8}  {'|d| (Å)':>10}"
              f"  {'Δx':>10}  {'Δy':>10}  {'Δz':>10}")
        for lid, dx, dy, dz, d in sorted_stats[:top_n]:
            print(f"  {lid:>5}  {local_to_lammps[lid]:>8}  {d:>10.6f}"
                  f"  {dx:>10.6f}  {dy:>10.6f}  {dz:>10.6f}")
        if args.verbose:
            for lid, dx, dy, dz, d in sorted_stats[top_n:]:
                print(f"  {lid:>5}  {local_to_lammps[lid]:>8}  {d:>10.6f}"
                      f"  {dx:>10.6f}  {dy:>10.6f}  {dz:>10.6f}")

    print()


if __name__ == '__main__':
    main()
