from celescope.tools import utils
from celescope.tools.star_mixin import Star_mixin, get_opts_star_mixin

class Star(Star_mixin):
    """
    ## Features
    - Align R2 reads to the reference genome with STAR.

    ## Output
    - `{sample}_Aligned.sortedByCoord.out.bam` BAM file contains Uniquely Mapped Reads.

    - `{sample}_SJ.out.tab` SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format.

    - `{sample}_Log.out` Main log with a lot of detailed information about the run. 
    This is most useful for troubleshooting and debugging.

    - `{sample}_Log.progress.out` Report job progress statistics, such as the number of processed reads, 
    % of mapped reads etc. It is updated in 1 minute intervals.

    - `{sample}_Log.Log.final.out` Summary mapping statistics after mapping job is complete, 
    very useful for quality control. The statistics are calculated for each read (single- or paired-end) and 
    then summed or averaged over all reads. Note that STAR counts a paired-end read as one read, 
    (unlike the samtools agstat/idxstats, which count each mate separately). 
    Most of the information is collected about the UNIQUE mappers 
    (unlike samtools agstat/idxstats which does not separate unique or multi-mappers). 
    Each splicing is counted in the numbers of splices, which would correspond to 
    summing the counts in SJ.out.tab. The mismatch/indel error rates are calculated on a per base basis, 
    i.e. as total number of mismatches/indels in all unique mappers divided by the total number of mapped bases.
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)


    
    @utils.add_log
    def run(self):
        super().run()



def star(args):
    with Star(args, display_title="Mapping") as runner:
        runner.run()

def get_opts_star(parser, sub_program):
    get_opts_star_mixin(parser, sub_program)
