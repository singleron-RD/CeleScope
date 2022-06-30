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
        self.count_fusion_df = None
        self.count_fusion_out = f'{self.out_prefix}_fusion_count.csv'

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
    
    def get_fusion_count_df(self):
        fusion_df = self.df_tsne.reset_index().filter(list(self.pos_dict.keys()) + ["barcode"])
        fusion_df = fusion_df.melt(id_vars=["barcode"],value_name="UMI",var_name="fusion",ignore_index=True)
        
        count = fusion_df[fusion_df["UMI"] > 0].groupby(["fusion"], as_index=False)["barcode"].count()
        self.count_fusion_df = count.assign(percent = lambda x: 100 * x.barcode / len(self.df_tsne))
        self.count_fusion_df.to_csv(self.count_fusion_out, index = None)
    
    def add_fusion_count_metrics(self):
        self.add_help_content('fusion_1, fusion_2...', 'number of positive cells(UMI > 0) of each fusion')
        for _, row in self.count_fusion_df.iterrows():
            self.add_metric(
                name=f'{row["fusion"]}',
                value=int(row['barcode']),
                total=int(len(self.df_tsne)),
            )

    def run(self):
        super().run()
        self.plot_fusion()
        self.get_fusion_count_df()
        self.add_fusion_count_metrics()