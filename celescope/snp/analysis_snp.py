import subprocess
import os
import sys
import warnings

import scanpy as sc
import pandas as pd
import pysam
import matplotlib
import matplotlib.pyplot as plt

from celescope.tools import utils
from celescope.tools.step import Step
from celescope.tools.step import s_common
from celescope.tools.target_metrics import get_genes
from celescope.__init__ import HELP_DICT
from celescope.snp.__init__ import PANEL


matplotlib.use("Agg")
warnings.filterwarnings("ignore")

AA_DICT = {
    "Gly": "G",
    "Ala": "A",
    "Val": "V",
    "Leu": "L",
    "Ile": "I",
    "Phe": "F",
    "Trp": "W",
    "Tyr": "Y",
    "Asp": "D",
    "Asn": "N",
    "Glu": "E",
    "Lys": "K",
    "Gln": "Q",
    "Met": "M",
    "Ser": "S",
    "Thr": "T",
    "Cys": "C",
    "Pro": "P",
    "His": "H",
    "Arg": "R",
}


def parse_variant_ann(variant_ann_file):
    """
    Args:
        variant_ann_file: variant annotation file from snpEff.

    Returns:
        gene_list, mRNA_list, protein_list
    """
    gene_list, mRNA_list, protein_list = [], [], []

    with open(variant_ann_file) as f:
        for line in f.readlines():
            if not line.startswith("#"):
                info = line.split("\t")[7]
                anns = info.split("|")
                gene = anns[3]
                gene_list.append(gene)

                tmp1, tmp2 = [], []
                for ann in anns:
                    if ann.startswith("c."):
                        exon_loc = anns[anns.index(ann) - 1].split("/")[0]
                        # WARNING_TRANSCRIPT_INCOMPLETE
                        if not exon_loc:
                            continue

                        exon = ann.strip("c.")
                        exon = f"exon{exon_loc}:{exon}"
                        if exon not in tmp1:
                            tmp1.append(exon)

                    if ann.startswith("p."):
                        protein = ann[2:]
                        for i in AA_DICT:
                            protein = protein.replace(i, AA_DICT[i])
                        if protein not in tmp2:
                            tmp2.append(protein)

                mRNA_list.append(",".join(tmp1))
                protein_list.append(",".join(tmp2))

    return (gene_list, mRNA_list, protein_list)


def parse_vcf_to_df(vcf_file, cols=("chrom", "pos", "alleles"), infos=("VID", "CID")):
    """
    Read cols and infos into pandas df
    """
    vcf = pysam.VariantFile(vcf_file)
    df = pd.DataFrame(columns=[col.capitalize() for col in cols] + infos)
    rec_dict = {}
    for rec in vcf.fetch():
        for col in cols:
            rec_dict[col.capitalize()] = getattr(rec, col)
            if col == "alleles":
                rec_dict["Alleles"] = "-".join(rec_dict["Alleles"])

        for info in infos:
            rec_dict[info] = rec.info[info]

        """
        rec_dict['GT'] = [s['GT'] for s in rec.samples.values()][0]
        rec_dict['GT'] = [str(item) for item in rec_dict['GT']]
        rec_dict['GT'] = '/'.join(rec_dict['GT'])
        """
        df_new = pd.DataFrame(rec_dict, index=[0])
        df = pd.concat([df, df_new])

    vcf.close()
    df.reset_index(drop=True, inplace=True)
    return df


def vcf_to_gt_csv(vcf_file, csv_file):
    vcf = pysam.VariantFile(vcf_file)

    samples = vcf.header.samples

    with open(csv_file, "w") as f:
        header = ["variant"] + list(samples)
        f.write(",".join(header) + "\n")

        for record in vcf:
            mutation_name = f"{record.chrom}_{record.pos}"
            genotypes = []

            for sample in samples:
                genotype = record.samples[sample]["GT"]
                g1, g2 = genotype

                if g1 is None:
                    genotype_str = "NA"
                else:
                    genotype_str = "/".join([str(g1), str(g2)])

                genotypes.append(genotype_str)

            line = [mutation_name] + genotypes
            f.write(",".join(line) + "\n")


class Analysis_snp(Step):
    """
    ## Features
    - Annotate variants with [snpEff](http://pcingola.github.io/SnpEff/).

    ## Output
    - `{sample}_gt.csv` Genotypes of variants of each cell. Rows are variants and columns are cells.
    - `{sample}_variant_ncell.csv` Number of cells with each genotype.
    - `{sample}_variant_table.csv` annotated with snpEff.

    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.vcf_file = args.vcf

        # parse
        self.genes = get_genes(args)
        self.database = args.database

        # data
        self.variant_table = None

        # out
        self.snpeff_outdir = f"{self.outdir}/snpEff/"
        self.snpeff_ann_vcf_file = f"{self.snpeff_outdir}/variants_ann.vcf"
        self.final_vcf_file = f"{self.out_prefix}_final.vcf"
        utils.check_mkdir(self.snpeff_outdir)
        self.plot_snp_dir = f"{self.outdir}/{self.sample}_plot_snp/"

        self.gt_file = f"{self.out_prefix}_gt.csv"
        self.ncell_file = f"{self.out_prefix}_variant_ncell.csv"
        self.variant_table_file = f"{self.out_prefix}_variant_table.csv"

    @utils.add_log
    def write_gt(self):
        vcf_to_gt_csv(self.final_vcf_file, self.gt_file)

    @utils.add_log
    def write_ncell(self):
        """
        parse gt_file to collect each genotype cell count into ncell_file
        """
        df = pd.read_csv(self.gt_file, index_col=0)
        df_ncell = df.apply(pd.Series.value_counts, axis=1).fillna(0).astype(int)
        df_ncell.to_csv(self.ncell_file, index=True)

    @utils.add_log
    def run_snpEff(self):
        args = f"-Xmx8g -v {self.database} {os.path.abspath(self.vcf_file)} "
        if self.args.snpeff_cache:
            args += f" -dataDir {self.args.snpeff_cache} "
        cmd = f"snpEff {args} > variants_ann.vcf "
        self.run_snpEff.logger.info(cmd)

        cwd = os.getcwd()
        os.chdir(self.snpeff_outdir)
        subprocess.check_call(cmd, shell=True)
        # change dir back to avoid can not find '09.analysis_snp/stat.txt' error
        os.chdir(cwd)

    @utils.add_log
    def keep_in_gene(self):
        """
        Output:
            self.final_vcf_file
        """
        gene_list, _, _ = parse_variant_ann(self.snpeff_ann_vcf_file)
        with pysam.VariantFile(self.snpeff_ann_vcf_file) as vcf_in:
            with pysam.VariantFile(
                self.final_vcf_file, "w", header=vcf_in.header
            ) as vcf_out:
                for i, record in enumerate(vcf_in.fetch()):
                    if gene_list[i] in self.genes:
                        vcf_out.write(record)

    def get_variant_table(self):
        """
        Returns:
            is_in_gene_list: if res[i] == True, line i is in gene_list
        """

        df_vcf = parse_vcf_to_df(self.final_vcf_file, infos=[])
        df_vcf["Gene"], df_vcf["mRNA"], df_vcf["Protein"] = parse_variant_ann(
            self.final_vcf_file
        )
        df_ncell = pd.read_csv(self.ncell_file)
        df_vcf = pd.concat([df_vcf, df_ncell], axis=1)

        cols = [
            "Chrom",
            "Pos",
            "Alleles",
            "Gene",
            "0/0",
            "0/1",
            "1/1",
            "mRNA",
            "Protein",
        ]
        cols = [col for col in cols if col in df_vcf.columns]
        df_vcf = df_vcf.loc[:, cols]
        is_in_gene_list = df_vcf.Gene.isin(self.genes)
        df_vcf = df_vcf[is_in_gene_list]

        self.variant_table = df_vcf
        self.variant_table.reset_index(drop=True, inplace=True)
        self.variant_table.to_csv(self.variant_table_file, index=False)

    def add_help(self):
        """
        <p> Chrom : chromosome name.</p>
        <p> Pos : the 1-based position of the variation on the given sequence..</p>
        <p> Alleles : REF(reference base or bases in the case of an indel) - ALT(alternative alleles).</p>
        <p> 0/0, 0/1, 1/1: number of cells with each genotype.</p>
        <p> Gene : gene symbol.</p>
        <p> mRNA :  A standard nomenclature is used in specifying the sequence changes.</p>
        <p> Protein :  A standard nomenclature is used in specifying the sequence changes.</p>
        """
        self.add_help_content(name="Chrom", content="Chromosome name")
        self.add_help_content(
            name="Pos",
            content="the 1-based position of the variation on the given sequence",
        )
        self.add_help_content(
            name="Alleles",
            content="REF(reference base or bases in the case of an indel) - ALT(alternative alleles)",
        )
        self.add_help_content(
            name="0/0, 0/1, 1/1", content="number of cells with each genotype"
        )
        self.add_help_content(name="Gene", content="gene symbol")
        self.add_help_content(
            name="mRNA",
            content="A standard nomenclature is used in specifying the sequence changes",
        )
        self.add_help_content(
            name="Protein",
            content="A standard nomenclature is used in specifying the sequence changes",
        )

    @utils.add_log
    def plot_snp(self):
        match_dict = utils.parse_match_dir(self.args.match_dir)
        if "h5ad" not in match_dict:
            return

        utils.check_mkdir(self.plot_snp_dir)
        df_gt = pd.read_csv(self.gt_file, keep_default_na=False, index_col=0)
        df_v = self.variant_table.copy()
        df_v["n_variants"] = df_v["0/1"] + df_v["1/1"]
        indices = df_v.nlargest(self.args.plot_top_n, "n_variants").index
        df_top = df_gt.iloc[indices,]
        df_top = df_top.transpose()
        variants = df_top.columns
        for c in variants:
            df_top[c] = df_top[c].astype("category")

        adata = sc.read_h5ad(match_dict["h5ad"])
        adata.obs = pd.concat([adata.obs, df_top], axis=1)
        pt_size = min(100, 120000 / len(adata.obs))
        gene_list, protein_list = df_v["Gene"], df_v["Protein"]
        for i, v in enumerate(variants):
            title = f"top{i+1}_{variants[i]}_{gene_list[indices[i]]}_{protein_list[indices[i]]}"
            file_name = f"{self.plot_snp_dir}/{title}.pdf"
            sc.pl.umap(
                adata,
                color=v,
                size=pt_size,
                palette={
                    "0/0": "dimgray",
                    "0/1": "orange",
                    "1/1": "red",
                    "NA": "lightgray",
                },
                title=title,
            )
            plt.savefig(file_name, dpi=300, bbox_inches="tight")

    def run(self):
        self.run_snpEff()
        self.keep_in_gene()
        self.write_gt()
        self.write_ncell()
        self.get_variant_table()
        self.add_help()
        try:
            self.plot_snp()
        except Exception as e:
            sys.stderr.write(f"{e}\n")
        table_dict = self.get_table_dict(
            title="Variant table", table_id="variant", df_table=self.variant_table
        )
        self.add_data(table_dict=table_dict)


@utils.add_log
def analysis_snp(args):
    with Analysis_snp(args, display_title="Analysis") as runner:
        runner.run()


def get_opts_analysis_snp(parser, sub_program):
    parser.add_argument("--gene_list", help=HELP_DICT["gene_list"])
    parser.add_argument(
        "--database",
        help="snpEff database. Common choices are GRCh38.mane.1.0.ensembl(human) and GRCm38.99(mouse)",
        default="GRCh38.mane.1.0.ensembl",
    )
    parser.add_argument("--snpeff_cache", help="snpEff cache directory.")
    parser.add_argument("--bed", help="custom bed file.")
    parser.add_argument("--panel", help=HELP_DICT["panel"], choices=list(PANEL))
    parser.add_argument(
        "--plot_top_n", type=int, help="plot UMAP of at most n variants ", default=20
    )
    if sub_program:
        s_common(parser)
        parser.add_argument("--match_dir", help=HELP_DICT["match_dir"], required=True)
        parser.add_argument("--vcf", help="vcf file.", required=True)
