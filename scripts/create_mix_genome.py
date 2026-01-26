import pysam
import argparse


def process_fasta(prefix1, fasta1, prefix2, fasta2):
    print("Processing FASTA files")
    out_fa = f"{prefix1}_{prefix2}.fasta"
    with open(out_fa, "w") as f_out:
        print(f"    - Processing species 1: {fasta1}")
        with pysam.FastxFile(fasta1) as fh1:
            for entry in fh1:
                f_out.write(f">{prefix1}_{entry.name}\n{entry.sequence}\n")
        print(f"    - Processing species 2: {fasta2}")
        with pysam.FastxFile(fasta2) as fh2:
            for entry in fh2:
                f_out.write(f">{prefix2}_{entry.name}\n{entry.sequence}\n")


def process_gtf(prefix1, gtf1, prefix2, gtf2):
    print("Processing GTF files...")
    out_gtf = f"{prefix1}_{prefix2}.gtf"
    with open(out_gtf, "w") as g_out:
        g_out.write(f"#mixed gtf of {prefix1} and {prefix2}\n")
        g_out.write(f"#{gtf1}\n")
        g_out.write(f"#{gtf2}\n")
        print(f"    - Processing species 1 GTF: {gtf1}")
        with open(gtf1, "r") as g1:
            for line in g1:
                if line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) >= 9:
                    cols[0] = f"{prefix1}_{cols[0]}"
                    g_out.write("\t".join(cols))

        print(f"    - Processing species 2 GTF: {gtf2}")
        with open(gtf2, "r") as g2:
            for line in g2:
                if line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) >= 9:
                    cols[0] = f"{prefix2}_{cols[0]}"
                    g_out.write("\t".join(cols))


def main():
    parser = argparse.ArgumentParser(
        description="Combine genomic references with species prefixes."
    )

    parser.add_argument("--prefix1", required=True, help="Prefix for species 1")
    parser.add_argument("--fasta1", required=True, help="Path to species 1 FASTA")
    parser.add_argument("--gtf1", required=True, help="Path to species 1 GTF")

    parser.add_argument("--prefix2", required=True, help="Prefix for species 2")
    parser.add_argument("--fasta2", required=True, help="Path to species 2 FASTA")
    parser.add_argument("--gtf2", required=True, help="Path to species 2 GTF")

    parser.add_argument("--out_fa", default="combined.fasta", help="Output FASTA name")

    args = parser.parse_args()

    process_fasta(args.prefix1, args.fasta1, args.prefix2, args.fasta2)
    process_gtf(args.prefix1, args.gtf1, args.prefix2, args.gtf2)

    print("Mix complete.")


if __name__ == "__main__":
    main()
