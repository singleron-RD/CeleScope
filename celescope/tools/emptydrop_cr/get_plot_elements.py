import collections
import math

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as pltoff

from celescope.tools.utils import add_log

BarcodeRankPlotSegment = collections.namedtuple('BarcodeRankPlotSegment', ['start', 'end', 'cell_density', 'legend'])

CHARTS_PLOTLY_MODEBAR_TRANSFORM_BUTTONS = [
    'zoom2d',
    'pan2d',
    'zoomIn2d',
    'zoomOut2d',
    'autoScale2d',
    # 'resetScale2d'  can't totally disable interaction, it seems-- keep reset option
]

CHARTS_PLOTLY_EXPORT_BUTTONS = [
    'toImage',
    'sendDataToCloud',
]

CHARTS_PLOTLY_FIXED_CONFIG = {
    'modeBarButtonsToRemove': CHARTS_PLOTLY_MODEBAR_TRANSFORM_BUTTONS+CHARTS_PLOTLY_EXPORT_BUTTONS,
    'displaylogo': False,
    'showLink': False
}

CHARTS_PLOTLY_MOVABLE_CONFIG = {
    'modeBarButtonsToRemove': CHARTS_PLOTLY_EXPORT_BUTTONS,
    'displaylogo': False,
    'showLink': False
}
BC_RANK_PLOT_LINE_WIDTH = 3
# Gradient scheme used in the barcode rank plot
BC_PLOT_COLORS = ['#dddddd', '#d1d8dc', '#c6d3dc', '#bacfdb', '#aecada', '#a3c5d9', '#97c0d9', '#8cbbd8', '#80b7d7',
                  '#74b2d7', '#6aadd6', '#66abd4', '#62a8d2', '#5ea5d1', '#59a2cf', '#559fce', '#519ccc', '#4d99ca',
                  '#4997c9', '#4594c7', '#4191c5', '#3d8dc4', '#3a8ac2', '#3787c0', '#3383be', '#3080bd', '#2c7cbb',
                  '#2979b9', '#2676b7', '#2272b6', '#1f6eb3', '#1d6ab0', '#1a65ac', '#1861a9', '#155ca6', '#1358a2',
                  '#10539f', '#0e4f9b', '#0b4a98', '#094695', '#09438f', '#0a4189', '#0c3f83', '#0d3d7c', '#0e3b76',
                  '#103970', '#11366a', '#123463', '#14325d', '#153057']

CHARTS = [
    {
        'layout': {
            'title': 'Barcode Rank',
            'width': 470,
            'height': 313,
            'margin': {'l': 60, 'r': 0, 't': 30, 'b': 40},
            'hovermode': 'closest',
            'xaxis': {
                'title': 'Barcodes',
                'type': 'log',
            },
            'yaxis': {
                'title': 'UMI counts',
                'type': 'log',
            },
        },
        'data': [],
        'config': CHARTS_PLOTLY_FIXED_CONFIG,
        'function': '_plot_barcode_rank',
        'description': 'The plot shows the count of filtered UMIs mapped to each barcode.  Barcodes can be determined to be cell-associated based on their UMI counts or by their expression profiles. Therefore some regions of the graph can contain both cell-associated and background-associated barcodes. The color of the graph represents the local density of barcodes that are cell-associated.',
        'name': 'barcode_rank',
    },
]


def BC_PLOT_CMAP(density):
    """
    Colormap utility fn to map a number to one of the colors in the gradient
    color scheme defined above
    Input
    - density : A real number in the range [0,1]
    """
    assert density >= 0.
    assert density <= 1.
    levels = len(BC_PLOT_COLORS)
    ind = min(levels - 1, int(math.floor(levels * density)))
    return BC_PLOT_COLORS[ind]


def get_plot_segment(start_index, end_index, sorted_bc, cell_barcodes, legend=False):
    """
    Helper function to build a plot segment.
    """
    assert end_index > start_index
    num_cells = sum([1 for i in range(start_index, end_index) if sorted_bc[i] in cell_barcodes])
    density = float(num_cells)/float(end_index-start_index)
    return BarcodeRankPlotSegment(start=start_index, end=end_index, cell_density=density, legend=legend)


def segment_log_plot_by_length(y_data, x_start, x_end):
    """
    Given the extends of the mixed region [x_start, x_end), compute
    the x-indices that would divide the plot into segments of a
    prescribed length (in pixel coordinates) with a minimum number
    of barcodes in each segment
    """

    if x_end <= x_start:
        return []

    SEGMENT_NORMALIZED_MAX_LEN = 0.02
    MIN_X_SPAN = 20

    log_max_x = np.log(len(y_data))
    log_max_y = np.log(max(y_data))

    this_segment_len = 0.0
    segment_idx = [x_start]

    for i in range(x_start, x_end):
        last_i = max(x_start, i-1)
        dx = (np.log(i) - np.log(last_i)) / log_max_x
        dy = (np.log(y_data[i]) - np.log(y_data[last_i])) / log_max_y
        this_segment_len += np.linalg.norm([dx, dy])
        if this_segment_len >= SEGMENT_NORMALIZED_MAX_LEN and i > (segment_idx[-1] + MIN_X_SPAN):
            segment_idx.append(i+1)
            this_segment_len = 0.0

    if segment_idx[-1] != x_end:
        segment_idx.append(x_end)

    return segment_idx


def convert_numpy_array_to_line_chart(array, ntype):
    array = np.sort(array)[::-1]

    rows = []
    previous_count = None
    for (index,), count in np.ndenumerate(array):
        if index == 0 or index == len(array)-1:
            rows.append([index, ntype(count)])
        elif previous_count != count:
            previous_index = rows[-1][0]
            if previous_index != index - 1:
                rows.append([index - 1, ntype(previous_count)])
            rows.append([index, ntype(count)])
        previous_count = count
    return rows


@add_log
def counter_barcode_rank_plot_data(count_data_path):
    """
    get cell density for each plot_segments
    :param count_data_path:
    :return: sorted_counts, plot_segments, cell_nums
    """
    count_data = pd.read_csv(count_data_path, index_col=0, sep='\t')
    cell_bc = np.array(count_data[count_data['mark'] == 'CB'].index)
    sorted_bc = np.array(count_data.index)
    sorted_counts = np.array(count_data['UMI'])
    cell_nums = len(cell_bc)
    total_bc = len(sorted_bc)
    # find the first barcode which is not a cell
    first_non_cell = total_bc
    for i, bc in enumerate(sorted_bc):
        if bc not in cell_bc:
            first_non_cell = i
            break

    # find the last barcode which is a cell
    last_cell = 0
    for i in reversed(range(total_bc)):
        if sorted_bc[i] in cell_bc:
            last_cell = i
            break

    ranges = [0, first_non_cell, last_cell+1, total_bc]
    plot_segments = []
    plot_segments.append(BarcodeRankPlotSegment(start=0, end=ranges[1], cell_density=1.0, legend=True))
    plot_segments.append(BarcodeRankPlotSegment(start=ranges[2], end=ranges[3], cell_density=0.0, legend=True))

    mixed_segments = segment_log_plot_by_length(sorted_counts, ranges[1], ranges[2])
    for i in range(len(mixed_segments) - 1):
        plot_segments.append(
            get_plot_segment(mixed_segments[i], mixed_segments[i + 1], sorted_bc, cell_bc, legend=False))

    return sorted_counts, plot_segments, cell_nums


def _plot_barcode_rank(chart, counts, num_cells):
    """ Generate a generic barcode rank plot """
    rows = convert_numpy_array_to_line_chart(counts, int)

    for _i, row in enumerate(rows):
        index, count = row[0], row[1]
        if index < num_cells:
            series_list = [chart['data'][0]]
        elif index == num_cells:
            # Connect the two lines
            series_list = [chart['data'][0], chart['data'][1]]
        else:
            series_list = [chart['data'][1]]

        for series in series_list:
            series['x'].append(index)
            series['y'].append(count)

    # Handle case where there is no background
    bg_series = chart['data'][1]
    if len(bg_series['x']) == 0:
        bg_series['x'].append(0)
        bg_series['y'].append(0)

    return chart


def build_plot_data_dict(plot_segment, counts):
    """
    Construct the data for a plot segment by appropriately slicing the
    counts
    Inputs:
    - plot_segment: BarcodeRankPlotSegment containing [start, end)
        of the segment, the cell density and legend visibility option
    - counts: Reverse sorted UMI counts for all barcodes.
    """

    start = max(0, plot_segment.start - 1)  # -1 for continuity between two charts
    end = plot_segment.end
    plot_rows = convert_numpy_array_to_line_chart(counts[start:end], int)
    name = 'Cells' if plot_segment.cell_density > 0 else 'Background'

    # Setup the tooltip
    if plot_segment.cell_density > 0.:
        n_barcodes = plot_segment.end - plot_segment.start
        n_cells = int(round(plot_segment.cell_density * n_barcodes))
        hover = "{:.0f}% Cells<br>({}/{})".format(100 * plot_segment.cell_density, n_cells, n_barcodes)
    else:
        hover = "Background"

    data_dict = {
        "x": [],
        "y": [],
        "name": name,
        "hoverinfo": "text",
        "text": hover,
        "type": "scattergl",
        "mode": "lines",
        "line": {
            "width": 3,
            "color": BC_PLOT_CMAP(plot_segment.cell_density),
        },
        "showlegend": plot_segment.legend,
    }
    offset = 1 + start  # it's a log-log plot, hence the 1
    for index, count in plot_rows:
        data_dict["x"].append(index + offset)
        data_dict["y"].append(count)

    # Handle case where the data is empty
    if len(data_dict["x"]) == 0:
        data_dict["x"].append(0)
        data_dict["y"].append(0)

    return data_dict


def _plot_counter_barcode_rank(chart, counts, plot_segments):
    """
    Generate the RNA counter barcode rank plot
    Inputs:
        - chart: chart element to populate data
        - counts: UMI counts reverse sorted
        - plot_segments: A list of BarcodeRankPlotSegments
    """

    for segment in plot_segments:
        chart["data"].append(build_plot_data_dict(segment, counts))

    return chart


@add_log
def get_plot_data(plot_segments, counts):
    plot_data = []
    for segment in plot_segments:
        plot_data.append(build_plot_data_dict(segment, counts))

    return plot_data


def plot_barcode_rank(count_file_path):
    sorted_counts, plot_segments, _cell_nums = counter_barcode_rank_plot_data(count_file_path)
    plot_data = get_plot_data(plot_segments, sorted_counts)

    plotly_data = [go.Scatter(x=dat['x'], y=dat['y'], name=dat['name'], mode=dat['mode'], showlegend=dat['showlegend'],
                              marker={'color': dat['line']['color']}, line=dat['line'], text=dat['text']) for dat in
                   plot_data]

    layout = go.Layout(width=470, height=313,
                       title={"text": "Barcode rank", "font": {"color": "black"}, "x": 0.5},
                       xaxis={"type": "log", "title": "Barcode", "titlefont": {"color": "black"},
                              "color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
                       yaxis={"type": "log", "title": "UMI counts", "titlefont": {"color": "black"},
                              "color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
                       margin=dict(l=50, r=0, t=30, b=30),
                       plot_bgcolor="#FFFFFF")

    config = dict({"displayModeBar": True,
                   "staticPlot": False,
                   "showAxisDragHandles": False,
                   "modeBarButtons": [["toImage", "resetScale2d"]],
                   "scrollZoom": False,
                   "displaylogo": False})

    fig = go.Figure(data=plotly_data, layout=layout)

    chart = pltoff.plot(fig, include_plotlyjs=True, output_type='div', config=config)

    return chart
