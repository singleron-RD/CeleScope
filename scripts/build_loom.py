#!/usr/bin/env python3

import argparse
import pandas as pd
import anndata
import scanpy as sc
import os


def build_loom_from_celescope(path, output):
    """
    Build a Loom file from celescope output files.

    Parameters
    ----------
    path : str
        Celescope output directory.
        Expected structure:
        path/
            ├── 01.starsolo/
            │     └── <prefix>_Solo.out/
            │           ├── Gene/raw/
            │           └── Velocyto/raw/
            └── outs/filtered/barcodes.tsv.gz

    output : str
        Output directory for loom file.

    Returns
    -------
    None
    """

    prefix = os.path.basename(os.path.abspath(path))

    starsolo_path = os.path.join(path, "01.starsolo", f"{prefix}_Solo.out")
    barcodes_used_path = os.path.join(path, "outs", "filtered", "barcodes.tsv.gz")

    print(f"[INFO] Processing sample: {prefix}")

    # === Read matrices ===
    print("[INFO] Reading matrices...")

    X = sc.read_mtx(os.path.join(starsolo_path, "Gene/raw/matrix.mtx")).X.transpose()
    spliced = sc.read_mtx(
        os.path.join(starsolo_path, "Velocyto/raw/spliced.mtx")
    ).X.transpose()
    unspliced = sc.read_mtx(
        os.path.join(starsolo_path, "Velocyto/raw/unspliced.mtx")
    ).X.transpose()
    ambiguous = sc.read_mtx(
        os.path.join(starsolo_path, "Velocyto/raw/ambiguous.mtx")
    ).X.transpose()

    # === Read metadata ===
    print("[INFO] Reading barcodes and features...")

    obs = pd.read_csv(
        os.path.join(starsolo_path, "Velocyto/raw/barcodes.tsv"),
        header=None,
        index_col=0,
    )

    var = pd.read_csv(
        os.path.join(starsolo_path, "Velocyto/raw/features.tsv"),
        sep="\t",
        names=("gene_ids", "feature_types"),
        index_col=1,
    )

    # === Create AnnData ===
    print("[INFO] Building AnnData object...")

    adata = anndata.AnnData(
        X=X,
        obs=obs,
        var=var,
        layers={"spliced": spliced, "unspliced": unspliced, "ambiguous": ambiguous},
    )

    adata.obs.index.names = ["CellID"]
    adata.var.index.names = ["Gene"]
    adata.var_names_make_unique()

    # === Filter cells ===
    print("[INFO] Filtering cells...")

    cells_used = pd.read_table(barcodes_used_path, header=None, index_col=0)

    adata = adata[list(cells_used.index), :]

    # === Add prefix to cell IDs ===
    adata.obs.index = [f"{prefix}:{cell}x" for cell in adata.obs.index]

    # === Write loom ===
    os.makedirs(output, exist_ok=True)
    loompath = os.path.join(output, f"{prefix}.loom")

    print("[INFO] Writing loom file...")
    adata.write_loom(loompath)

    print(f"[DONE] Loom file saved at: {loompath}")


def main():
    parser = argparse.ArgumentParser(
        description="Build loom file from CeleScope STARsolo Velocyto output"
    )

    parser.add_argument(
        "-i", "--input", required=True, help="Celescope output directory"
    )

    parser.add_argument("-o", "--output", required=True, help="Output directory")

    args = parser.parse_args()

    build_loom_from_celescope(args.input, args.output)


if __name__ == "__main__":
    main()
