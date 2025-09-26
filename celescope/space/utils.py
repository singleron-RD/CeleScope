import pandas as pd
from pathlib import Path


class Spatial:
    """
    take a Abrowser spatial as input and output a visium like spatial directory
    """

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.tissue_positions_list = self.input_dir / "tissue_positions_list.csv"
        self.scalefactors = self.input_dir / "scalefactors_json.json"
        self.tissue_loweres = self.input_dir / "tissue_lowres.png"
        self.tissue_hires = self.input_dir / "tissue_hires.png"

    def get_in_tissue_barcodes(self) -> list[str]:
        # 读取 CSV，没有表头
        positions = pd.read_csv(self.tissue_positions_list, header=None)

        # 设置列名
        positions.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
        ]

        barcodes_in_tissue = positions.loc[
            positions["in_tissue"] == 1, "barcode"
        ].tolist()
        return barcodes_in_tissue
