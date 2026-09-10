---
layout: default
title: apply_event.py
parent: Scripts
nav_order: 1
---

# apply_event.py
{: .no_toc }

Applique un événement ARTn à une configuration LAMMPS pour produire les configurations selle et finale.
{: .fs-5 .fw-300 }

---

## Description

Lit trois fichiers d'entrée — une configuration LAMMPS initiale, un fichier d'événement ARTn identifiant les atomes du cluster, et un fichier XYZ à trois frames (initial, selle, final) — calcule les déplacements atomiques et écrit deux nouveaux fichiers LAMMPS avec les conditions aux limites périodiques (PBC) appliquées par défaut.

**Dépendances :** bibliothèque standard Python 3 uniquement.

---

## Usage

```bash
./apply_event.py init.data event_file event.xyz [options]

# avec un nom de base personnalisé
./apply_event.py conf.data event_001 event_001.xyz -o run_001

# sans PBC, afficher les 20 plus grands déplacements
./apply_event.py conf.data event_001 event_001.xyz --no-pbc --top 20
```

---

## Arguments

| Argument | Statut | Description |
|:---|:---|:---|
| `init.data` | Requis | Configuration LAMMPS initiale (`.data`) |
| `event_file` | Requis | Fichier d'événement ARTn — IDs des atomes du cluster |
| `event.xyz` | Requis | Fichier XYZ à 3 frames : initial, selle, final |
| `-o, --output` | Optionnel | Nom de base des fichiers de sortie |
| `--saddle-out` | Optionnel | Chemin explicite pour la configuration selle |
| `--final-out` | Optionnel | Chemin explicite pour la configuration finale |
| `--no-pbc` | Optionnel | Désactiver l'encapsulation PBC |
| `-v, --verbose` | Optionnel | Afficher les déplacements par atome |
| `--top N` | Optionnel | Afficher les N plus grands déplacements (défaut : 10) |

---

## Fichiers de sortie

- `<base>_saddle.data` — configuration à l'état selle
- `<base>_final.data` — configuration à l'état final
