---
layout: default
title: neb_conf2lammps.py
nav_order: 3
---

# neb_conf2lammps.py
{: .no_toc }

Convertit des configurations NEB (Nudged Elastic Band) en format de données LAMMPS.
{: .fs-5 .fw-300 }

---

## Description

Fusionne la géométrie de boîte et les identifiants d'atomes de la configuration NEB initiale avec les positions atomiques de la configuration NEB finale, et produit un fichier `.data` compatible LAMMPS. Vérifie que les deux fichiers décrivent le même nombre d'atomes.

**Dépendances :** bibliothèque standard Python 3 uniquement.

---

## Usage

```bash
./neb_conf2lammps.py neb_ini.conf neb_fin.conf

# avec nom de sortie explicite
./neb_conf2lammps.py neb_ini.conf neb_fin.conf -o final_state.data
```

---

## Arguments

| Argument | Statut | Description |
|:---|:---|:---|
| `ini_file` | Requis | Configuration NEB initiale (format LAMMPS) |
| `fin_file` | Requis | Configuration NEB finale avec les positions cibles |
| `-o, --output` | Optionnel | Chemin de sortie (défaut : `fin_file` avec extension `.data`) |
