from collections import defaultdict

import plotly
import plotly.express as px
import plotly.graph_objects as go

from celescope.tools import utils


PLOTLY_CONFIG = {
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
        "t": 30, }
}


class Plotly_plot():

    def __init__(self, df):
        self._df = df
        self._fig = None

    def plotly_plot(self):
        return plotly.offline.plot(
            self._fig,
            include_plotlyjs=True,
            output_type='div',
            config=PLOTLY_CONFIG
        )

    def get_plotly_div(self):

        return self.plotly_plot()


class Tsne_plot(Plotly_plot):

    def __init__(self, df_tsne, feature_name, discrete=True):
        super().__init__(df_tsne)

        self.feature_name = feature_name
        self.discrete = discrete
        title_feature_name = feature_name[0].upper() + feature_name[1:]
        self.title = f"t-SNE plot Colored by {title_feature_name}"

        self._layout = {}
        self._dot_size = 4
        self._df['size'] = self._dot_size
        self._df['barcode_index'] = list(range(1, len(self._df) + 1))
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

    def set_color_scale(self, color_scale):
        self.scatter_config['color_continuous_scale'] = color_scale

    @utils.add_log
    def get_plotly_div(self):

        if self.discrete:
            self.discrete_tsne_plot()
        else:
            self.continuous_tsne_plot()
        self.update_fig()

        return self.plotly_plot()

    @utils.add_log
    def discrete_tsne_plot(self):

        sum_df = self._df.groupby([self.feature_name]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res_dict = defaultdict(int)
        res_list = []
        for cluster in sorted(self._df[self.feature_name].unique()):
            name = f"{cluster}({percent_df[cluster]}%)"
            res_dict[cluster] = name
            res_list.append(name)

        self._df[self.feature_name] = self._df[self.feature_name].map(res_dict)

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
            title={"text": self.title, "x": 0.5, "y": 0.95, "font": {"size": 15}},
            plot_bgcolor='#FFFFFF',
            hovermode="closest"
        )



class Tsne_dropdown_plot(Plotly_plot):

    def __init__(self, df_tsne,name, feature_name_list, dropdown=True):
        super().__init__(df_tsne)
        
        self.name = name
        self.feature_name_list = feature_name_list
        self.dropdown = dropdown
        self.title = f"t-SNE plot Colored by {self.name}"

        self._layout = {}
        
        self._buttons=[]
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
        
        
    def get_plotly_div(self):
        
        if self.dropdown:
            self._fig = go.Figure()
            _num = len(self.feature_name_list)
            i=0
            for feature_name in self.feature_name_list:
                self._fig.add_trace(go.Scatter(x=self._df[self._str_coord1],y=self._df[self._str_coord2],mode='markers',
                                               name = feature_name[0].upper() + feature_name[1:],
                                               showlegend=False,
                                               marker=go.scatter.Marker(opacity=0.9,size=4,color=self._df[feature_name],
                                                                        colorscale=[[0,'LightGrey'],[1,'Red']],
                                                                        colorbar=go.scatter.marker.ColorBar(title=feature_name,
                                                                                                            titlefont={'size':11})),
                                               textposition='top center'
                                                    ))
                visible_list=[False]*_num
                visible_list[i]=True
                i += 1
                button = dict(args=[{"visible":visible_list},],
                  label=feature_name,
                  method="update",
                 )
                self._buttons.append(button)
        else:
            pass
        self.update_fig()
        
        return self.plotly_plot()
        
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
            updatemenus=[{'buttons':self._buttons,
                          'direction':'down',
                          'pad':{'r':8,'t':8},
                          'showactive':True,
                          'x':0.9,
                          'xanchor':'left',
                          'y':1.2,
                          'yanchor':'top'}],
            annotations=[{'text': 'TAG type:',
                          'showarrow': False,
                          'x': 0.9,
                          'xref': 'paper',
                          'xanchor': 'right',
                          'y': 1.15,
                          'yref': 'paper',
                          'font': {'size': 11},
                          'align': 'left'}],
            title={"text": self.title, "x": 0.5, "y": 0.95, "font": {"size": 15}},
            plot_bgcolor='#FFFFFF',
            hovermode="closest"
        )
        
        self._fig.update_yaxes(
            title_text=self._str_coord2,
            **self.axes_config
        )
        
        self._fig.update_yaxes(
            title_text=self._str_coord2,
            **self.axes_config
        )


class Bar_plot(Plotly_plot):

    def __init__(self, df_bar):
        super().__init__(df_bar)
        self.set_fig()
    
    def set_fig(self):
        self._fig = px.bar(
                    x=[str(i) for i in list(self._df.head(10).ClonotypeID)],
                    y=self._df.head(10).proportion.tolist(),
                    labels={'x': 'Clonotype ID', 'y': 'Fraction of Cells'},
                    width=700,
                    height=500
                    )
        self._fig.update_traces(marker_color='rgb(158,202,225)', marker_line_color='rgb(8,48,107)',
                   marker_line_width=1.5, opacity=0.6)
        self._fig.update_layout(title_text='Top 10 Clonotype Frequencies',
                   title={"x": 0.5, "y": 0.9, "font": {"size": 20, "family": "San Serif"}},
                   plot_bgcolor='#FFFFFF')

                   
class Pie_plot(Plotly_plot):

    def __init__(self, df_region):
        super().__init__(df_region)
        self.set_fig()

    def set_fig(self):
        layout = {
            "height": 300,
            "width": 400,
            "margin": {
                "l": 50,
                "r": 35,
                "b": 10,
                "t": 10,
            }
        }
        self._fig = px.pie(
            self._df,
            names='regions',
            values='values',
        )
        self._fig.update_traces(textposition='none')
        self._fig.update_layout(layout)

class Line_plot(Plotly_plot):

    def __init__(self, df_line,title=None,x_title=None,y_title=None,y_range=None,section=True,):
        super().__init__(df_line)
        self.df_line = df_line
        self.title = title
        self.x_title = x_title
        self.y_title = y_title
        self.y_range = y_range
        self.section = section

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
        
        if self.section:
            self.yaxes_config = {
                'showgrid': True,
                'gridcolor': '#F5F5F5',
                'linecolor':'black',
                'showline': True, 
                'ticks': None,
                'rangemode':'tozero',
            }
        else:
            self.yaxes_config = {
                'showgrid': True,
                'gridcolor': '#F5F5F5',
                'linecolor':'black',
                'showline': True, 
                'ticks': None,
                'range':self.y_range,
            }

        self.line_config = {
            'data_frame': df_line,
            'title': self.title,
            'x': self.x_title, 
            'y': self.y_title,
        }

        self.line_plot()
        self.update_fig()

    @utils.add_log
    def line_plot(self):
        self._fig = px.line(
            **self.line_config,
        )  

    def update_fig(self):
        self._fig.update_xaxes(
            **self.xaxes_config
        )

        self._fig.update_yaxes(
            **self.yaxes_config
        )
        
        self._fig.update_layout(
            LAYOUT,
            title={"x":0.5, "y":0.95, "font":{"size":15}},
            yaxis_zeroline = True,
            showlegend = False,
            plot_bgcolor = '#FFFFFF',
            hovermode = "closest"
        )
