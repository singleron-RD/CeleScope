import math
import numpy as np
from PIL import Image

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
    "displaylogo": False,
}

COLORS = px.colors.qualitative.Plotly + px.colors.qualitative.Alphabet

LAYOUT = {
    "height": 313,
    "width": 400,
    "margin": {
        "l": 45,
        "r": 35,
        "b": 30,
        "t": 30,
    },
}


class Plotly_plot:
    def __init__(self, df):
        self._df = df
        self._fig = None

    def plotly_plot(self):
        return plotly.offline.plot(
            self._fig, include_plotlyjs=False, output_type="div", config=PLOTLY_CONFIG
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
        self._df["size"] = self._dot_size
        self._df["barcode_index"] = list(range(1, len(self._df) + 1))
        self._str_coord1 = "tSNE_1"
        self._str_coord2 = "tSNE_2"
        self.axes_config = {
            "showgrid": True,
            "gridcolor": "#F5F5F5",
            "showline": False,
            "ticks": None,
            "zeroline": True,
            "zerolinecolor": "black",
            "zerolinewidth": 0.7,
        }

        self.scatter_config = {
            "data_frame": df_tsne,
            "title": self.title,
            "x": self._str_coord1,
            "y": self._str_coord2,
            "size_max": self._dot_size,
            "hover_data": {
                self._str_coord1: False,
                self._str_coord2: False,
                self.feature_name: True,
                "barcode_index": True,
                "size": False,
            },
            "size": "size",
            "opacity": 0.9,
            "color": self.feature_name,
            "color_discrete_sequence": COLORS,
            "color_continuous_scale": px.colors.sequential.Jet,
        }

    def set_color_scale(self, color_scale):
        self.scatter_config["color_continuous_scale"] = color_scale

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
        res_dict = dict()
        res_list = []
        for cluster in sorted(self._df[self.feature_name].unique()):
            name = f"{cluster}({percent_df[cluster]}%)"
            res_dict[cluster] = name
            res_list.append(name)

        self._df[self.feature_name] = self._df[self.feature_name].map(res_dict)

        self._fig = px.scatter(
            **self.scatter_config, category_orders={self.feature_name: res_list}
        )

    @utils.add_log
    def continuous_tsne_plot(self):
        self._fig = px.scatter(
            **self.scatter_config,
        )

    def update_fig(self):
        self._fig.update_xaxes(title_text=self._str_coord1, **self.axes_config)

        self._fig.update_yaxes(title_text=self._str_coord2, **self.axes_config)

        self._fig.update_layout(
            self._layout,
            title={"text": self.title, "x": 0.5, "y": 0.95, "font": {"size": 15}},
            plot_bgcolor="#FFFFFF",
            hovermode="closest",
        )


class Tsne_dropdown_plot(Plotly_plot):
    def __init__(self, df_tsne, name, feature_name_list):
        super().__init__(df_tsne)

        self.name = name
        self.feature_name_list = feature_name_list
        self.title = f"t-SNE plot Colored by {self.name}"
        self.x_pos = 0.1
        self._layout = {}

        self._buttons = []
        self.x_label, self.y_label = "tSNE_1", "tSNE_2"

        self.x = df_tsne[self.x_label]
        self.y = df_tsne[self.y_label]
        self.axes_config = {
            "showgrid": True,
            "gridcolor": "#F5F5F5",
            "showline": False,
            "ticks": None,
            "zeroline": True,
            "zerolinecolor": "black",
            "zerolinewidth": 0.7,
        }

        x_ = math.ceil(max(abs(self.x)))
        y_ = math.ceil(max(abs(self.y)))
        self.x_range = [-x_, x_]
        self.y_range = [-y_, y_]

    def get_plotly_div(self):
        self.continuous_tsne_plot()
        self.update_fig()
        return self.plotly_plot()

    def continuous_tsne_plot(self):
        self._fig = go.Figure()

        # remove the cluster fraction "cluster 1(10%)"
        self._df["hover_text"] = "cluster " + self._df["cluster"].apply(
            lambda x: x.split("(")[0]
        )

        # Coordinates
        for i, feature in enumerate(self.feature_name_list):
            self._fig.add_trace(
                go.Scattergl(
                    x=self.x,
                    y=self.y,
                    mode="markers",
                    marker=dict(
                        size=4,
                        color=self._df[feature],
                        colorbar=dict(title="UMI"),
                        colorscale=px.colors.sequential.Jet,
                        opacity=0.9,
                    ),
                    name=feature,
                    hovertext=self._df["hover_text"],  # Text to display when hovering
                    hoverinfo="text",  # Display only the hover text on hover
                    visible=True if i == 0 else False,
                )
            )

    def update_fig(self):
        self._fig.update_xaxes(
            range=self.x_range, title_text=self.x_label, **self.axes_config
        )

        self._fig.update_yaxes(
            range=self.y_range, title_text=self.y_label, **self.axes_config
        )

        # Create dropdown buttons to toggle traces
        buttons = []
        for i, feature in enumerate(self.feature_name_list):
            visibility = [False] * len(self.feature_name_list)
            visibility[i] = True
            buttons.append(
                dict(
                    label=feature,
                    method="update",
                    args=[
                        {"visible": visibility},
                    ],
                )
            )

        self._fig.update_layout(
            self._layout,
            updatemenus=[
                {
                    "buttons": buttons,
                    "direction": "down",
                    "pad": {"r": 8, "t": 8},
                    "showactive": True,
                    "x": 1.1,
                    "xanchor": "right",
                    "y": 1.2,
                    "yanchor": "top",
                }
            ],
            annotations=[
                {
                    "text": "TAG type:",
                    "showarrow": False,
                    "x": self.x_pos,
                    "xref": "paper",
                    "xanchor": "right",
                    "y": 1.15,
                    "yref": "paper",
                    "font": {"size": 11},
                    "align": "left",
                }
            ],
            title={"text": self.title, "x": 0.5, "y": 0.95, "font": {"size": 15}},
            plot_bgcolor="#FFFFFF",
            hovermode="closest",
        )


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
            },
        }
        self._fig = px.pie(
            self._df,
            names="regions",
            values="values",
        )
        self._fig.update_traces(textposition="none")
        self._fig.update_layout(layout)


class Line_plot(Plotly_plot):
    def __init__(
        self,
        df_line,
        title=None,
        x_title=None,
        y_title=None,
        y_range=None,
        section=True,
    ):
        super().__init__(df_line)
        self.df_line = df_line
        self.title = title
        self.x_title = x_title
        self.y_title = y_title
        self.y_range = y_range
        self.section = section

        self.xaxes_config = {
            "showgrid": True,
            "gridcolor": "#F5F5F5",
            "linecolor": "black",
            "showline": True,
            "ticks": None,
            "tickmode": "linear",
            "tick0": 0,
            "dtick": 0.5,
        }

        if self.section:
            self.yaxes_config = {
                "showgrid": True,
                "gridcolor": "#F5F5F5",
                "linecolor": "black",
                "showline": True,
                "ticks": None,
                "rangemode": "tozero",
            }
        else:
            self.yaxes_config = {
                "showgrid": True,
                "gridcolor": "#F5F5F5",
                "linecolor": "black",
                "showline": True,
                "ticks": None,
                "range": self.y_range,
            }

        self.line_config = {
            "data_frame": df_line,
            "title": self.title,
            "x": self.x_title,
            "y": self.y_title,
        }

        self.line_plot()
        self.update_fig()

    @utils.add_log
    def line_plot(self):
        self._fig = px.line(
            **self.line_config,
        )

    def update_fig(self):
        self._fig.update_xaxes(**self.xaxes_config)

        self._fig.update_yaxes(**self.yaxes_config)

        self._fig.update_layout(
            LAYOUT,
            title={"x": 0.5, "y": 0.95, "font": {"size": 15}},
            yaxis_zeroline=True,
            showlegend=False,
            plot_bgcolor="#FFFFFF",
            hovermode="closest",
        )


class Bar_plot(Plotly_plot):
    def __init__(self, df_bar):
        super().__init__(df_bar)
        self.set_fig()

    def set_fig(self):
        self._fig = px.bar(
            x=[str(i) for i in list(self._df.head(10).ClonotypeID)],
            y=self._df.head(10).proportion.tolist(),
            labels={"x": "Clonotype ID", "y": "Fraction of Cells"},
            width=700,
            height=500,
        )
        self._fig.update_traces(
            marker_color="rgb(158,202,225)",
            marker_line_color="rgb(8,48,107)",
            marker_line_width=1.5,
            opacity=0.6,
        )
        self._fig.update_layout(
            title_text="Top 10 Clonotype Frequencies",
            title={"x": 0.5, "y": 0.9, "font": {"size": 20, "family": "San Serif"}},
            plot_bgcolor="#FFFFFF",
        )


class Conversion_plot(Plotly_plot):
    def __init__(self, df_bar):
        super().__init__(df_bar)
        np.random.seed(1)
        self.values, bins = np.histogram(df_bar["cell_pct"], bins=200)
        self.bin_centers = (bins[:-1] + bins[1:]) / 2
        self.set_fig()

    def set_fig(self):
        self._fig = px.bar(x=self.bin_centers, y=self.values)
        self._fig.update_layout(yaxis_type="log")
        self._fig.update_xaxes(title="cell percent")
        self._fig.update_yaxes(title="count")
        self._fig.update_layout(plot_bgcolor="rgba(0,0,0,0)")


class Substitution_plot(Plotly_plot):
    def __init__(self, df_bar):
        super().__init__(df_bar)
        self.set_fig()

    def set_fig(self):
        num4sample = 0
        colors4sample = {}
        num4x = 0

        self._fig = go.Figure()

        for sample in self._df["sample"].unique():
            legend_show = True
            colors4sample[sample] = COLORS[num4sample]
            num4sample += 1
            flag_x = "x" + str(num4x + 1)
            df_plot = self._df[self._df["sample"] == sample]
            num4x += 1

            self._fig.add_trace(
                go.Bar(
                    name=sample + "+",
                    x=df_plot["sample"],
                    y=df_plot["+"],
                    legendgroup=sample,
                    marker_color=colors4sample[sample],
                    marker_line_color="#FFFFFF",
                    showlegend=legend_show,
                    xaxis=flag_x,
                )
            )
            self._fig.add_trace(
                go.Bar(
                    name=sample + "-",
                    x=df_plot["sample"],
                    y=df_plot["-"],
                    legendgroup=sample,
                    showlegend=legend_show,
                    marker_color=colors4sample[sample],
                    marker_line_color="#FFFFFF",
                    opacity=0.3,
                    xaxis=flag_x,
                )
            )

        self._fig.update_layout(barmode="stack")

        per = 1 / (num4x + 1)
        gap4bar = per / len(self._df["sample"].unique())
        num4x = 0
        for _ in self._df["sample"].unique():
            if num4x == 0:
                flag_x = "xaxis"
            else:
                flag_x = "xaxis" + str(num4x + 1)
            anchor_x = "x" + str(num4x + 1)
            num4x += 1
            self._fig["layout"][flag_x] = dict(
                domain=[per * num4x, per * (num4x + 1) - gap4bar], anchor=anchor_x
            )

        self._fig.update_layout(plot_bgcolor="#FFFFFF")
        self._fig.update_xaxes(
            showgrid=False,
            linecolor="black",
            showline=True,
            ticks="outside",
            showticklabels=True,
        )
        self._fig.update_yaxes(
            showgrid=False, linecolor="black", showline=True, ticks="outside"
        )
        # width_num = 400 * (len(self._df['sample'].unique()) * len(self._df['sample'].unique())) / (5*12)  # set bar width
        self._fig.update_layout(height=500)  # , width=width_num)
        self._fig.update_layout(legend=dict(orientation="h"))
        self._fig.update_layout(
            legend=dict(
                yanchor="top",
                y=1.3,
                xanchor="left",
                x=0.07,
                valign="top",
            )
        )
        self._fig.update_layout(yaxis_title="Rates of nucleotide substitution (%)")


class Violin_plot(Plotly_plot):
    def __init__(self, df_vln, feature_name, color="#1f77b4", label=True):
        super().__init__(df_vln)
        self.feature_name = feature_name
        self.color = color
        self.label = label

        if self.label:
            self.y = "Fraction of labeled UMI"
        else:
            self.y = "Fraction of new RNA"
        self.set_fig()

    def set_fig(self):
        self._fig = go.Figure()
        self._fig.add_trace(
            go.Violin(
                y=self._df,
                box_visible=True,
                meanline_visible=True,
                line_color="black",
                fillcolor=self.color,
                opacity=0.6,
                x0=self.feature_name,
            )
        )
        self._fig.update_layout(yaxis_zeroline=True, showlegend=False)
        self._fig.update_layout(plot_bgcolor="#FFFFFF")
        self._fig.update_xaxes(
            showgrid=False, linecolor="black", showline=True, ticks=None
        )
        self._fig.update_yaxes(
            showgrid=False,
            linecolor="black",
            showline=True,
            ticks="outside",
            title_text=self.y,
            rangemode="tozero",
        )


class StaticPlot:
    def __init__(self, img_path):
        img = Image.open(img_path)
        if img.mode != "RGB":
            img = img.convert("RGB")
        fig = go.Figure()
        width, height = img.size
        # 添加图片
        fig.add_layout_image(
            dict(
                source=img,
                xref="x",
                yref="y",
                x=0,
                y=height,
                sizex=width,
                sizey=height,
                sizing="stretch",
                opacity=1,
                layer="below",
            )
        )

        # 隐藏坐标轴
        fig.update_xaxes(visible=False, range=[0, width], fixedrange=False)
        fig.update_yaxes(
            visible=False,
            range=[0, height],
            scaleanchor="x",
            autorange="reversed",
            fixedrange=False,
        )

        # 设置 Figure 尺寸，并去掉边距
        fig.update_layout(
            autosize=True,
            dragmode="pan",
            margin=dict(l=0, r=0, t=0, b=0),  # 左右上下边距都设为0
            paper_bgcolor="white",  # 整个图纸背景
            plot_bgcolor="white",
        )
        self._fig = fig

    def get_div(self):
        return plotly.offline.plot(
            self._fig,
            include_plotlyjs=False,
            output_type="div",
            config=dict(scrollZoom=True),
        )
