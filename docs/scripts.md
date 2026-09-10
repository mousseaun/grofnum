---
layout: default
title: Scripts
nav_order: 2
has_children: true
---

# Scripts

Scripts utilitaires pour le post-traitement de simulations LAMMPS, ARTn et KMC.

| Script | Rôle | Dépendances |
|:---|:---|:---|
| [`apply_event.py`](apply_event) | Applique un événement ARTn à une configuration LAMMPS | stdlib |
| [`neb_conf2lammps.py`](neb_conf2lammps) | Convertit des configurations NEB en format LAMMPS | stdlib |
| [`show_basin_connectivity.py`](show_basin_connectivity) | Affiche le graphe de connectivité des bassins pykmc | pandas |
| [`view_trajectory.py`](view_trajectory) | Visualiseur 3-D interactif de trajectoires KMC | OVITO, PySide6 |
| [`xyz2lammps.py`](xyz2lammps) | Convertit XYZ étendu en fichier de données LAMMPS | stdlib |
| [`xyz_distance.py`](xyz_distance) | Calcule le RMSD entre deux configurations atomiques | stdlib |
| [`nettoyage.sh`](nettoyage) | Supprime les fichiers de sortie LAMMPS/pykmc | — |
