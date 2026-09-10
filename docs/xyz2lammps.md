---
layout: default
title: xyz2lammps.py
parent: Scripts
nav_order: 5
---

# xyz2lammps.py
{: .no_toc }

Convertit des fichiers XYZ étendus en fichiers de données LAMMPS (`atom_style atomic`).
{: .fs-5 .fw-300 }

---

## Description

Parse les en-têtes XYZ étendus pour extraire la géométrie de la boîte de simulation (paramètres tricliniques inclus), mappe les symboles élémentaires en types LAMMPS avec les masses atomiques standards, et supporte les fichiers de trajectoire multi-frames via `--frame`.

**Dépendances :** bibliothèque standard Python 3 uniquement.

---

## Usage

```bash
# conversion simple (sortie : input.data)
./xyz2lammps.py structure.xyz

# avec symboles et nom de sortie
./xyz2lammps.py structure.xyz -o al2o3.data --symbols Al,O

# frame 5 d'une trajectoire
./xyz2lammps.py traj.xyz --frame 5

# lister toutes les frames disponibles
./xyz2lammps.py traj.xyz --list-frames
```

---

## Arguments

| Argument | Statut | Description |
|:---|:---|:---|
| `input` | Requis | Fichier XYZ étendu (structure unique ou trajectoire) |
| `-o, --output` | Optionnel | Fichier de sortie (défaut : `input` avec extension `.data`) |
| `--frame N` | Optionnel | Frame à convertir (défaut : 0) |
| `--list-frames` | Optionnel | Lister toutes les frames et quitter |
| `--symbols` | Optionnel | Symboles élémentaires séparés par virgule, ex. `Si` ou `Al,O` |
