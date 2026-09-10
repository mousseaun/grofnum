---
layout: default
title: Accueil
nav_order: 1
description: "Scripts utilitaires pour LAMMPS, ARTn et KMC — mousseaun/grofnum"
permalink: /
---

# grofnum
{: .no_toc }

Scripts utilitaires pour le post-traitement de simulations de dynamique moléculaire (LAMMPS), d'événements ARTn, et de trajectoires KMC.
{: .fs-5 .fw-300 }

---

## Scripts disponibles

| Script | Rôle | Dépendances |
|:---|:---|:---|
| [`apply_event.py`](apply_event) | Applique un événement ARTn à une configuration LAMMPS | stdlib |
| [`neb_conf2lammps.py`](neb_conf2lammps) | Convertit des configurations NEB en format LAMMPS | stdlib |
| [`show_basin_connectivity.py`](show_basin_connectivity) | Affiche le graphe de connectivité des bassins pykmc | pandas |
| [`view_trajectory.py`](view_trajectory) | Visualiseur 3-D interactif de trajectoires KMC | OVITO, PySide6 |
| [`xyz2lammps.py`](xyz2lammps) | Convertit XYZ étendu en fichier de données LAMMPS | stdlib |
| [`xyz_distance.py`](xyz_distance) | Calcule le RMSD entre deux configurations atomiques | stdlib |
| [`nettoyage.sh`](nettoyage) | Supprime les fichiers de sortie LAMMPS/pykmc | — |

---

## Dépendances

La majorité des scripts n'utilisent que la bibliothèque standard Python 3 (`argparse`, `math`, `re`, `sys`, `os`). Les exceptions sont :

- **`show_basin_connectivity.py`** — nécessite `pandas`
- **`view_trajectory.py`** — nécessite `ovito` et `PySide6`

## Installation

```bash
git clone https://github.com/mousseaun/grofnum.git
cd grofnum
# rendre les scripts exécutables (optionnel)
chmod +x *.py nettoyage.sh
```
