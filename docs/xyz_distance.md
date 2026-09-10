---
layout: default
title: xyz_distance.py
nav_order: 7
---

# xyz_distance.py
{: .no_toc }

Calcule la distance structurale (RMSD) entre deux configurations atomiques.
{: .fs-5 .fw-300 }

---

## Description

Accepte indifféremment des fichiers XYZ étendus et des fichiers de données LAMMPS (format détecté automatiquement à partir de l'extension, surchargeable). Les atomes sont appariés par champ `id` s'il est présent, sinon par ordre de ligne. Rapporte le RMSD en Å, le RMS par composante (dx, dy, dz), le déplacement maximal, et la différence d'énergie si elle est disponible dans les deux fichiers.

**Dépendances :** bibliothèque standard Python 3 uniquement.

---

## Usage

```bash
# comparer deux fichiers LAMMPS
./xyz_distance.py initial.data final.data

# comparer des frames spécifiques d'une trajectoire XYZ
./xyz_distance.py traj.xyz ref.xyz --frame-a 10 --frame-b 0

# mode verbeux : afficher les déplacements par atome, top 5
./xyz_distance.py A.data B.data -v --top 5

# désactiver PBC (cluster isolé)
./xyz_distance.py A.xyz B.xyz --no-pbc
```

---

## Arguments

| Argument | Statut | Description |
|:---|:---|:---|
| `fichier_A` | Requis | Première configuration — XYZ ou LAMMPS data |
| `fichier_B` | Requis | Deuxième configuration — même format ou différent |
| `--frame-a N` / `--frame-b N` | Optionnel | Index de frame à lire dans chaque fichier XYZ |
| `--format-a` / `--format-b` | Optionnel | Forcer le format : `xyz` ou `lammps` |
| `--no-pbc` | Optionnel | Désactiver les conditions aux limites périodiques |
| `-v, --verbose` | Optionnel | Afficher les vecteurs de déplacement par atome |
| `--top N` | Optionnel | Afficher les N atomes les plus déplacés |

---

## Sortie

- RMSD (Å)
- Déplacement maximum et atome concerné
- RMS par composante (dx, dy, dz)
- ΔE (si présent dans les deux fichiers)
