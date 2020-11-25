import os
import glob 
import argparse
import logging


class Debug():
    def __init__(self):
        self.outdir = 'debug'
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

    def opts(self):
        readme = f'debug'
        parser = argparse.ArgumentParser(readme)
        parser.add_argument('--dir', required=True)
        parser.add_argument('--n_sub', default=1000000)
        parser.add_argument('--genomeDir', default='/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92')
        args = parser.parse_args()
        self.dir = args.dir
        self.n_sub = args.n_sub
        self.genomeDir = args.genomeDir
        self.thread = 2

    def run_subsample(self):
        self.all_reads = glob.glob(f'{self.dir}/02.cutadapt/*clean_2.fq.gz')[0]
        self.sub_reads = f'{self.outdir}/sub.fq.gz'
        cmd = (
            f'less {self.all_reads} '
            f'|head -n {self.n_sub * 4} | gzip -c > {self.sub_reads} '
        )
        with open('subsample.sh', 'wt') as f:
            f.write(cmd)

    def run_STAR(self):
        cmd = (
            f'STAR --runThreadN {self.thread} \\\n'
            f'--genomeDir {self.genomeDir} \\\n'
            f'--readFilesIn {self.sub_reads} \\\n'
            f'--readFilesCommand zcat --outFilterMultimapNmax 1 \\\n'
            f'--outFileNamePrefix {self.outdir}/sub_ \\\n'
            f'--outSAMtype BAM SortedByCoordinate --outFilterMatchNmin 0 \\\n'
            f'--outReadsUnmapped Fastx \\\n'
        )
        with open('STAR.sh', 'wt') as f:
            f.write(cmd)
        self.unmapped = f'{self.outdir}/sub_Unmapped.out.mate1'

    def fastqc(self):
        cmd = (
            f'fastqc {self.unmapped} \\\n'
        )
        with open('fastqc.sh', 'wt') as f:
            f.write(cmd)


    def run(self):
        self.opts()
        self.run_subsample()
        self.run_STAR()
        self.fastqc()


if __name__ == '__main__':
    de = Debug()
    de.run()

    