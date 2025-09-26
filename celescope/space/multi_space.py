from collections import defaultdict
from celescope.tools import utils
from celescope.space.__init__ import __ASSAY__
from celescope.tools.multi import Multi, get_read


class Multi_space(Multi):
    @staticmethod
    @utils.add_log
    def parse_mapfile(mapfile, default_val):
        fq_dict = defaultdict(lambda: defaultdict(list))
        col4_dict = {}
        col5_dict = {}
        with open(mapfile) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                line_split = line.split()
                library_id, library_path, sample_name = line_split[:3]
                if len(line_split) >= 4:
                    col4 = line_split[3]
                else:
                    col4 = default_val
                fq1 = ",".join(get_read(library_id, library_path=library_path, read=1))

                fq_dict[sample_name]["fq1"].append(fq1)
                fq_dict[sample_name]["col4"].append(col4)
                col4_dict[sample_name] = col4
                if len(line_split) == 5:
                    col5_dict[sample_name] = line_split[4]

        for sample_name in fq_dict:
            fq_dict[sample_name]["fq1_str"] = ",".join(fq_dict[sample_name]["fq1"])

        if not fq_dict:
            raise Exception("empty mapfile!")
        return fq_dict, col4_dict, col5_dict

    def starsolo(self, sample):
        step = "starsolo"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} --spatial {self.col4_dict[sample]} '
        )
        self.process_cmd(
            cmd,
            step,
            sample,
            m=int(self.args.limitBAMsortRAM / 1e9),
            x=self.args.thread,
        )


def main():
    multi = Multi_space(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
