import argparse
import pysam


def main():

    parser = argparse.ArgumentParser('Add barcode and UMI tag to celescope BAM')
    parser.add_argument('--bam', help='input bam', required=True)
    parser.add_argument('--outdir', help='output directory', default='./')
    args = parser.parse_args()
    add_tag(args.bam, args.outdir)


def add_tag(bam, outdir):

    bam_out = f"{outdir}/add_tag.bam"
    handle_in = pysam.AlignmentFile(bam, "rb")
    header = handle_in.header
    handle_out = pysam.AlignmentFile(bam_out, "wb", header=header)
    for read in handle_in:
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        read.set_tag(tag='CB', value=barcode, value_type='Z')
        read.set_tag(tag='UB', value=umi, value_type='Z')
        handle_out.write(read)
    handle_in.close()
    handle_out.close()


if __name__ == '__main__':
    main()