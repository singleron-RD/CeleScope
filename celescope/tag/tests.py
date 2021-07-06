import subprocess
import unittest


class Tests(unittest.TestCase):
    def setUp(self):
        pass

    def test_seurat_hashtag(self):
        umi_tag = "/SGRNJ03/randd/RD20062501_yiqidemo/20210312/SMK-5C-BX-TAG/04.count_tag/SMK-5C-BX-TAG_umi_tag.tsv"
        matrix_10X = "/SGRNJ03/randd/RD20062501_yiqidemo/20210309/SMK-5C-BX-ZL/05.count/SMK-5C-BX-ZL_matrix_10X/"
        outdir = './'
        sample = 'test1'
        cmd = (
            'Rscript /SGRNJ/Database/script/pipe/develop/dev_CeleScope/celescope/tag/seurat_hashtag.R '
            f'--outdir {outdir} '
            f'--sample {sample} '
            f'--umi_tag {umi_tag} '
            f'--matrix_10X {matrix_10X} '
        )
        print(cmd)
        subprocess.check_call(cmd, shell=True)
