---
layout: default
title: nettoyage.sh
nav_order: 8
---

# nettoyage.sh
{: .no_toc }

Supprime les fichiers de sortie LAMMPS et pykmc du répertoire courant.
{: .fs-5 .fw-300 }

---

## Description

Script shell sans argument qui supprime définitivement les fichiers temporaires et de résultats générés par LAMMPS et pykmc dans le répertoire courant.

---

## Usage

```bash
# depuis votre répertoire de simulation
bash nettoyage.sh

# ou après chmod +x
./nettoyage.sh
```

---

## Fichiers supprimés

| Motif | Description |
|:---|:---|
| `artn.out*` | Fichiers de sortie ARTn |
| `lammps.log*` | Journaux LAMMPS |
| `pykmc.*` | Fichiers pykmc |
| `trajkmc.xyz` | Trajectoire KMC au format XYZ |
| `*.pickle` | Fichiers sérialisés Python (basin_connectivity, etc.) |

{: .warning }
Suppression définitive sans confirmation. Vérifiez le répertoire courant avant d'exécuter le script.
