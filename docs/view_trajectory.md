---
layout: default
title: view_trajectory.py
parent: Scripts
nav_order: 4
---

# view_trajectory.py
{: .no_toc }

Visualiseur 3-D interactif de trajectoires KMC via OVITO.
{: .fs-5 .fw-300 }

---

## Description

Charge un fichier XYZ multi-frames et le rend dans un viewport OVITO avec une interface Qt. Permet la navigation frame par frame ou la lecture automatique. Inclut une analyse Common Neighbor Analysis (CNA) pour identifier les défauts structuraux, avec masquage optionnel des atomes en structure FCC pour isoler les environnements non-FCC.

Il présente un exemple d'automatisation de la visualisation avec Ovito-python.

**Dépendances :** `ovito`, `PySide6`

---

## Usage

```bash
# fichier par défaut (trajkmc.xyz dans le répertoire courant)
python view_trajectory.py

# chemin explicite
python view_trajectory.py simulations/run_07/trajkmc.xyz

# afficher la documentation intégrée
python view_trajectory.py --info
```

---

## Arguments

| Argument | Statut | Description |
|:---|:---|:---|
| `fichier.xyz` | Optionnel | Chemin du fichier de trajectoire (défaut : `trajkmc.xyz`) |
| `--info` | Optionnel | Afficher la documentation et quitter |

---

## Interface

L'interface propose :

- Un curseur de navigation entre les frames
- Des boutons lecture / pause
- Le masquage des atomes FCC (pour isoler les défauts)
- L'affichage du nombre d'atomes et du nombre total de frames
