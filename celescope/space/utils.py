import shutil
import h5py
import numpy as np
import pandas as pd
from pathlib import Path
import scanpy as sc
import json


class Spatial:
    """
    take a Abrowser spatial as input and output a visium like spatial directory
    """

    def __init__(self, input_dir):
        self.input_dir = Path(input_dir)
        self.tissue_positions_list = self.input_dir / "tissue_positions_list.csv"
        self.scalefactors = self.input_dir / "scalefactors_json.json"
        self.tissue_loweres = self.input_dir / "tissue_lowres.png"
        self.tissue_hires = self.input_dir / "tissue_hires.png"
        self.parquet = self.input_dir / "tissue_positions.parquet"

        self.positions = pd.read_csv(self.tissue_positions_list, header=None)
        self.positions.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]

    def get_in_tissue_barcodes(self) -> list[str]:
        barcodes_in_tissue = self.positions.loc[
            self.positions["in_tissue"] == 1, "barcode"
        ].tolist()
        return barcodes_in_tissue

    def output_spatial(self, outdir):
        outdir = Path(outdir)
        shutil.copytree(self.input_dir, outdir, dirs_exist_ok=True)
        parquet = outdir / "tissue_positions.parquet"
        self.positions.to_parquet(parquet, index=False)

        # add "microns_per_pixel"，"regist_target_img_scalef" to scalefactors
        scalefactors = outdir / "scalefactors_json.json"
        with open(scalefactors) as f:
            data = json.load(f)
        data["microns_per_pixel"] = 1.0
        data["regist_target_img_scalef"] = 1.0
        with open(scalefactors, "w") as f:
            json.dump(data, f, indent=4)


def convert_10x_h5(mtx_dir, outfile, library_id="library0"):
    """
    Convert a 10X MTX folder to a Visium-compatible HDF5 file.
    Existing file will be overwritten.

    Args:
        mtx_dir (str or Path): Path to 10X-formatted mtx directory.
        outfile (str or Path): Output HDF5 file path.
        library_id (str): Library ID for Visium metadata.
    """
    mtx_dir = Path(mtx_dir)
    outfile = Path(outfile)
    if not outfile.suffix == ".h5":
        outfile = outfile.with_suffix(".h5")

    # 如果文件已存在，直接覆盖
    if outfile.exists():
        outfile.unlink()

    # 读取 MTX 文件
    adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)

    n_cells = adata.n_obs
    n_genes = adata.n_vars

    # 默认 feature_type 和 genome
    feature_types = (
        adata.var["feature_types"]
        if "feature_types" in adata.var.columns
        else ["Gene Expression"] * n_genes
    )
    genomes = (
        adata.var["genome"] if "genome" in adata.var.columns else ["unknown"] * n_genes
    )
    gene_ids = (
        adata.var["gene_ids"] if "gene_ids" in adata.var.columns else adata.var_names
    )

    # 创建 HDF5 文件
    with h5py.File(outfile, "w") as f:
        grp = f.create_group("matrix")

        # barcodes
        grp.create_dataset(
            "barcodes",
            data=np.array(adata.obs_names, dtype=f"|S{max(map(len, adata.obs_names))}"),
        )

        # counts
        grp.create_dataset("data", data=adata.X.data.astype(np.int32))
        grp.create_dataset("indices", data=adata.X.indices.astype(np.int32))
        grp.create_dataset("indptr", data=adata.X.indptr.astype(np.int32))
        grp.create_dataset("shape", data=np.array([n_genes, n_cells], dtype=np.int64))

        # features
        ftrs = grp.create_group("features")
        ftrs.create_dataset(
            "id", data=np.array(gene_ids, dtype=f"|S{max(map(len, gene_ids))}")
        )
        ftrs.create_dataset(
            "name",
            data=np.array(adata.var_names, dtype=f"|S{max(map(len, adata.var_names))}"),
        )
        ftrs.create_dataset(
            "feature_type",
            data=np.array(feature_types, dtype=f"|S{max(map(len, feature_types))}"),
        )
        ftrs.create_dataset(
            "genome", data=np.array(genomes, dtype=f"|S{max(map(len, genomes))}")
        )

        # Visium-specific metadata
        all_tag_keys = np.array(
            [b"feature_type", b"genome", b"id", b"name", b"library_ids"]
        )
        ftrs.create_dataset("_all_tag_keys", data=all_tag_keys)

        library_ids_dataset = np.array([library_id.encode()] * n_genes)
        ftrs.create_dataset("library_ids", data=library_ids_dataset)

        # HDF5 attributes required by scanpy.read_visium
        f.attrs["library_ids"] = np.array([library_id.encode()])
        f.attrs["chemistry_description"] = "Spatial3"
        f.attrs["software_version"] = "CustomPythonScript"
