import os
import unittest

import numpy as np
import pandas as pd

from .otsu import array2hist, makePlot, threshold_otsu


class test_capture(unittest.TestCase):
    def setUp(self):
        pass

    def test_otsu(self):
        count_files = [
            '/SGRNJ02/RandD4/virus_panel/20210124_4/S1225_EBV_Skin_Auto_SDF_NEB/07.count_virus/S1225_EBV_Skin_Auto_SDF_NEB_virus_UMI_count.tsv',
            '/SGRNJ02/RandD4/virus_panel/20210124/virus_test3_R_A_Beads_Manual_KZ/04.count_capture_virus/virus_test3_R_A_Beads_Manual_KZ_virus_UMI_count.tsv',
            '/SGRNJ02/RandD4/virus_panel/20210124/virus_test3_R_A_3Mins_Manual_KZ/04.count_capture_virus/virus_test3_R_A_3Mins_Manual_KZ_virus_UMI_count.tsv'
        ]
        count_file = '/SGRNJ02/RandD4/virus_panel/20210124_4/S1225_EBV_Skin_Auto_SDF_NEB/07.count_virus/S1225_EBV_Skin_Auto_SDF_NEB_virus_UMI_count.tsv'
        for count_file in count_files:
            df = pd.read_csv(count_file, sep='\t')
            array = np.log10(df["UMI"])
            hist = array2hist(array)
            print(hist)
            thresh = threshold_otsu(hist)
            fname = os.path.basename(count_file) + '.png'
            makePlot(hist, thresh, fname)
