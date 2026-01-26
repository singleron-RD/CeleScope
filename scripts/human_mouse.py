"""
Script to analyze human-mouse mixed single-cell RNA-seq data and identify doublets.
Generates a summary table with counts of human cells, mouse cells, doublets, and ambients.
Split the matrix into human and mouse based on gene names.
"""

import sys
from celescope.tools.matrix import CountMatrix
import pandas as pd
import argparse
import glob
import os
import functools


def run(mtx_path, doublet_threshold: float) -> dict:
    """
    Split human and mouse cells based on gene names and identify doublets.
    """
    mtx = CountMatrix.from_matrix_dir(mtx_path)
    sample = mtx_path.split("/")[-3]
    m = mtx.get_matrix()
    f = mtx.get_features()

    human_genes = []
    mouse_genes = []
    for i, gn in enumerate(f.gene_name):
        if gn == gn.upper():
            human_genes.append(i)
        else:
            mouse_genes.append(i)
    m_csr = m.tocsr()
    umi_h = m_csr[human_genes, :].sum(axis=0)
    umi_m = m_csr[mouse_genes, :].sum(axis=0)
    df_m = pd.DataFrame(umi_m.T, columns=["mouse"])
    df_h = pd.DataFrame(umi_h.T, columns=["human"])
    df = pd.merge(df_m, df_h, left_index=True, right_index=True)
    df["umi_sum"] = df.sum(axis=1)
    df["human_percent"] = df["human"] / df["umi_sum"] * 100
    df["mouse_percent"] = df["mouse"] / df["umi_sum"] * 100
    df["identity"] = "doublet"
    df.loc[df["human_percent"] < doublet_threshold, "identity"] = "mouse"
    df.loc[df["mouse_percent"] < doublet_threshold, "identity"] = "human"
    df.loc[df["identity"] == "human", "ambient_percent"] = df[
        df["identity"] == "human"
    ].mouse_percent
    df.loc[df["identity"] == "mouse", "ambient_percent"] = df[
        df["identity"] == "mouse"
    ].human_percent

    n_cell = df.shape[0]
    n_doublet = df[df["identity"] == "doublet"].shape[0]
    doublet_percent = round(n_doublet / n_cell * 100, 2)
    dict = {
        "n_cell": n_cell,
        "n_human_cell": df[df["identity"] == "human"].shape[0],
        "n_mouse_cell": df[df["identity"] == "mouse"].shape[0],
        "doublet_umi_percent_threshold": doublet_threshold,
        "n_doublet": n_doublet,
        "doublet_cell_percent": doublet_percent,
        "median_ambient_umi_percent": round(df["ambient_percent"].median(), 2),
        "mean_ambient_umi_percent": round(df["ambient_percent"].mean(), 2),
    }

    # split matrix
    human_barcodes = df[df["identity"] == "human"].index.tolist()
    mouse_barcodes = df[df["identity"] == "mouse"].index.tolist()
    out_dir = sample + ".out"
    mtx.slice_matrix(human_barcodes).to_matrix_dir(
        os.path.join(out_dir, "human_filtered")
    )
    mtx.slice_matrix(mouse_barcodes).to_matrix_dir(
        os.path.join(out_dir, "mouse_filtered")
    )

    return dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--celescope_dir",
        type=str,
        help="comma separated celescope output dirs",
        required=True,
    )
    parser.add_argument(
        "--doublet_threshold",
        type=float,
        default=25.0,
        help="threshold to define doublet cells",
    )
    args = parser.parse_args()

    cur_dir = os.getcwd()
    dfs = []
    matrix_path = set()
    for d in args.celescope_dir.split(","):
        paths = glob.glob(os.path.join(d, "*/*/filtered"))
        matrix_path.update(paths)
    if not matrix_path:
        sys.exit(
            "No filtered matrix found in the provided celescope directories. The filtered matrix should be located at <celescope_dir>/<sample>/outs/filtered."
        )
    print("Analyzing matrices from paths:", matrix_path)

    for path in matrix_path:
        sample = path.split("/")[-3]
        dict = run(path, args.doublet_threshold)
        dfs.append(pd.DataFrame.from_dict(dict, orient="index", columns=[sample]))

    df = functools.reduce(
        lambda left, right: pd.merge(left, right, left_index=True, right_index=True),
        dfs,
    )
    df = df.T
    df = df.astype(
        {
            "n_cell": int,
            "n_human_cell": int,
            "n_mouse_cell": int,
            "n_doublet": int,
        }
    )
    os.chdir(cur_dir)
    current_time = pd.Timestamp.now().strftime("%Y%m%dT%H%M%S")
    outprefix = f"human_mouse_summary-{current_time}"
    df.to_csv(f"{outprefix}.tsv", sep="\t")


if __name__ == "__main__":
    main()
