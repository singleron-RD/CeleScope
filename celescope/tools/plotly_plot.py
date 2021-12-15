from collections import defaultdict

import pandas as pd
import plotly
import plotly.express as px

import celescope.tools.utils as utils


PLOTLY_CONFIG =  {
    "displayModeBar": True, 
    "staticPlot": False, 
    "showAxisDragHandles": False, 
    "modeBarButtons": [["toImage", "resetScale2d"]], 
    "scrollZoom": False,
    "displaylogo": False
}

COLORS = px.colors.qualitative.Plotly + px.colors.qualitative.Alphabet

LAYOUT = {
        "height": 313,
        "width": 400,
        "margin": {
            "l": 45,
            "r": 35,
            "b": 30,
            "t": 30,}
}


class Tsne_plot():

    def __init__(self, df_tsne, feature_name, discrete=True):
        self.df_tsne = df_tsne
        self.feature_name = feature_name
        self.discrete = discrete
        title_feature_name = feature_name[0].upper() + feature_name[1:]
        self.title = f"t-SNE plot Colored by {feature_name}"
        
        self._layout = {}
        self._dot_size = 4
        self.df_tsne['size'] = self._dot_size
        self.df_tsne['barcode_index'] = list(range(1, len(self.df_tsne) + 1))
        self._str_coord1 = "tSNE_1"
        self._str_coord2 = "tSNE_2"
        self.axes_config = {
            'showgrid': True,
            'gridcolor': '#F5F5F5',
            'showline': False, 
            'ticks': None,
            'zeroline': True,
            'zerolinecolor': 'black',
            'zerolinewidth': 0.7,
        }

        self.scatter_config = {
            'data_frame': df_tsne,
            'title': self.title,
            'x': self._str_coord1, 
            'y': self._str_coord2,
            'size_max': self._dot_size,
            'hover_data': {
                self._str_coord1: False,
                self._str_coord2: False,
                self.feature_name: True,
                'barcode_index': True,
                'size': False,
            },
            'size': 'size',
            'opacity': 0.9,
            'color': self.feature_name,
            'color_discrete_sequence': COLORS,
            'color_continuous_scale': px.colors.sequential.Jet,
        }

        self._fig = None

        if discrete:
            self.discrete_tsne_plot()
        else:
            self.continuous_tsne_plot()
        self.update_fig()

        self.plotly_div = plotly.offline.plot(self._fig, include_plotlyjs=True, output_type='div', config=PLOTLY_CONFIG)

    @utils.add_log
    def discrete_tsne_plot(self):
        
        df_tsne = self.df_tsne
        feature_name = self.feature_name

        sum_df = df_tsne.groupby([feature_name]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res_dict = defaultdict(int)
        res_list = []
        for cluster in sorted(df_tsne[feature_name].unique()):
            name = f"{cluster}({percent_df[cluster]}%)"
            res_dict[cluster]= name
            res_list.append(name)

        df_tsne[self.feature_name] = df_tsne[self.feature_name].map(res_dict)

        self._fig = px.scatter(
            **self.scatter_config,
            category_orders={self.feature_name: res_list}
        )                               

    @utils.add_log
    def continuous_tsne_plot(self):

        self._fig = px.scatter(
            **self.scatter_config,
        )   

    def update_fig(self):
        self._fig.update_xaxes(
            title_text=self._str_coord1,
            **self.axes_config
        )

        self._fig.update_yaxes(
            title_text=self._str_coord2,
            **self.axes_config
        )
        
        self._fig.update_layout(
            self._layout,
            title={ "text":self.title, "x":0.5, "y":0.95, "font":{"size":15} },
            plot_bgcolor = '#FFFFFF',
            hovermode="closest"
        )

class Line_plot():
    def __init__(self, df_line,index: int):
        self.df_line = df_line
        self.index = index
        self.title = ['Sequencing Saturation','Median Genes per Cell']
        
        self._str_coord1 = "Reads Fraction"
        self._str_coord2 = ["Sequencing Saturation(%)","Median Genes per Cell"]

        self.xaxes_config = {
            'showgrid': True,
            'gridcolor': '#F5F5F5',
            'linecolor':'black',
            'showline': True, 
            'ticks': None,
            'tickmode':'linear',
            'tick0':0,
            'dtick':0.5,
        }

        self.yaxes_config = [{
            'showgrid': True,
            'gridcolor': '#F5F5F5',
            'linecolor':'black',
            'showline': True, 
            'ticks': None,
            'range':[0,100]},
            {
            'showgrid': True,
            'gridcolor': '#F5F5F5',
            'linecolor':'black',
            'showline': True, 
            'ticks': None,
            'rangemode':'tozero',
        }]

        self.line_config = {
            'data_frame': df_line,
            'title': self.title[index-1],
            'x': self._str_coord1, 
            'y': self._str_coord2[index-1],
        }

        self.line_plot()
        self.update_fig()

        self.plotly_div = plotly.offline.plot(self._fig, include_plotlyjs=True, output_type='div', config=PLOTLY_CONFIG)

    @utils.add_log
    def line_plot(self):
        self._fig = px.line(
            **self.line_config,
        )  

    def update_fig(self):
        index = self.index
        self._fig.update_xaxes(
            title_text=self._str_coord1,
            **self.xaxes_config
        )

        self._fig.update_yaxes(
            title_text=self._str_coord2[index-1],
            **self.yaxes_config[index-1]
        )
        
        self._fig.update_layout(
            LAYOUT,
            title={"x":0.5, "y":0.95, "font":{"size":15}},
            yaxis_zeroline = True,
            showlegend = False,
            plot_bgcolor = '#FFFFFF',
            hovermode = "closest"
        )
