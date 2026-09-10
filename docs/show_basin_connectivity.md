---
layout: default
title: show_basin_connectivity.py
parent: Scripts
nav_order: 3
---

# show_basin_connectivity.py
{: .no_toc }

Affiche et filtre le DataFrame de connectivité des bassins issu d'un fichier pickle pykmc.
{: .fs-5 .fw-300 }

---

## Description

Lit un fichier pickle pandas contenant le DataFrame `basin_connectivity` produit par pykmc. Supporte le filtrage par état source et état cible, ainsi que l'export CSV. Gère automatiquement les incompatibilités de versions pandas.

**Dépendances :** `pandas`

---

## Usage

```bash
# fichier par défaut (basin_connectivity_1.pickle)
python3 show_basin_connectivity.py

# fichier spécifique + filtre par état + export
python3 show_basin_connectivity.py run_42.pickle \
  --state 5 --to-state 3 --csv transitions_5to3.csv
```

---

## Arguments

| Argument | Statut | Description |
|:---|:---|:---|
| `pickle_file` | Optionnel | Chemin du fichier pickle (défaut : `basin_connectivity_1.pickle`) |
| `--state N` | Optionnel | Filtrer les lignes où `state` == N |
| `--to-state N` | Optionnel | Filtrer les lignes où `state_connexion` == N |
| `--csv FILE` | Optionnel | Exporter le tableau filtré en CSV |
