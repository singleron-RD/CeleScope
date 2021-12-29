import plotnine as p9

from celescope.tools.capture.analysis import Analysis, get_opts_analysis
from celescope.fusion.count_fusion import Count_fusion
from celescope.fusion.mkref import Mkref_fusion

def analysis_fusion(args):

    with Analysis_fusion(args) as runner:
        runner.run()


def get_opts_analysis_fusion(parser, sub_program):
    parser.add_argument('--fusion_genomeDir', help='Fusion genome directory.', required=True)
    get_opts_analysis(parser, sub_program)


class Analysis_fusion(Analysis):
    def __init__(self, args, display_title='Analysis'):
        super().__init__(args, display_title)

        fusion_pos_file =  Mkref_fusion.parse_genomeDir(args.fusion_genomeDir)['fusion_pos']
        self.pos_dict = Count_fusion.read_pos_file(fusion_pos_file)

        self.p9_theme = {    
            'axis_line_x': p9.element_line(size=2, colour="black"),
            'axis_line_y':p9.element_line(size=2, colour="black"),
            'panel_grid_major':p9.element_blank(),
            'panel_grid_minor':p9.element_blank(),
            'panel_border':p9.element_blank(),
            'panel_background':p9.element_blank(),
            'axis_text_x':p9.element_text(colour="black"),
            'axis_text_y':p9.element_text(colour="black"),
        }

    def plot_fusion(self):
        """
        plot fusion count
        """
        
        p9.theme_set(p9.theme_void())
        for ref in self.pos_dict:
            if ref in self.df_tsne.columns:
                out_plot_file = f'{self.out_prefix}_{ref}_fusion.pdf'
                plot = p9.ggplot(self.df_tsne, p9.aes(x="tSNE_1", y="tSNE_2", color=ref)) + \
                    p9.geom_point(size=0.2) + \
                    p9.theme_bw() + \
                    p9.scale_color_gradient(low="lightgrey",high="blue")
                plot.save(out_plot_file)

    def run(self):
        super().run()
        self.plot_fusion()





