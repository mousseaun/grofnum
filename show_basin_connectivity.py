#!/usr/bin/env python3
"""Display the content of a basin_connectivity pickle file (pandas DataFrame)."""

import argparse
import sys

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pickle_file", nargs="?", default="basin_connectivity_1.pickle",
                         help="Path to the basin_connectivity pickle file (default: %(default)s)")
    parser.add_argument("--state", type=int, default=None,
                         help="Only show rows where 'state' equals this value")
    parser.add_argument("--to-state", type=int, default=None,
                         help="Only show rows where 'state_connexion' equals this value")
    parser.add_argument("--csv", default=None,
                         help="Also export the table to this CSV file")
    args = parser.parse_args()

    try:
        df = pd.read_pickle(args.pickle_file)
    except NotImplementedError as exc:
        sys.exit(
            f"Failed to unpickle '{args.pickle_file}' with pandas {pd.__version__}: {exc}\n"
            "This usually means the file was written with a newer pandas that uses a "
            "different string dtype for column labels (e.g. pandas >= 3.0).\n"
            "Fix: upgrade pandas in this environment, e.g.\n"
            "    pip install -U pandas\n"
            "or read it once from an environment with a matching pandas version and "
            "re-save it (e.g. df.to_csv(...) or df.to_pickle(..., protocol=4))."
        )

    if args.state is not None:
        df = df[df["state"] == args.state]
    if args.to_state is not None:
        df = df[df["state_connexion"] == args.to_state]

    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)

    print(f"File: {args.pickle_file}")
    print(f"Shape: {df.shape[0]} rows x {df.shape[1]} columns")
    print(f"Columns: {list(df.columns)}")
    print()
    print(df.to_string(index=True))

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"\nSaved to {args.csv}")


if __name__ == "__main__":
    sys.exit(main())
