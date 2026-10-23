import ast
from copy import deepcopy
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_rgba
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import pytest

from tdpy.plotting import apply_plot_theme, configure_plot_style, plot_background_colors, plot_style


@pytest.fixture
def plot_functions(tmp_path):
    path = Path(__file__).with_name('plotter.py')
    tree = ast.parse(path.read_text())
    names = {'minor_tick_format', 'outcome_dist_hists', 'is_within_limits', 'scatter_plot'}
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]

    def late_style(axis):
        axis.set_facecolor('gray')
        axis.figure.set_facecolor('gray')

    namespace = {
        'np': np, 'pd': pd, 'math': math, 'plt': plt, 'deepcopy': deepcopy,
        'Normalize': Normalize, 'ScalarMappable': ScalarMappable, 'Line2D': Line2D, 'ticker': ticker,
        'plot_style': plot_style, 'configure_plot_style': configure_plot_style, 'apply_plot_theme': apply_plot_theme,
        'IsBurst': True, 'DATA_DIR': tmp_path, 'real_outliers': set(), 'adjust_text': None,
        'smplotlib': None, 'retro_noir': late_style,
        'Freq_1': pd.Series([20.0, 30.0]), 'L_EU_1': pd.Series([0.001, 0.002]),
        'Freq_2': pd.Series([120.0, 150.0], index=[2, 3]), 'L_EU_2': pd.Series([0.001, 0.002], index=[2, 3]),
        'NenuFreq': np.array([15.0, 85.0]), 'NenuNoise': np.array([0.001, 0.002]),
        'uGMRT': pd.DataFrame({'Frequencies': [[100.0, 200.0]] * 4, 'RMS Noise': [[0.001, 0.002]] * 4}),
        'MWA': pd.DataFrame({'Frequencies': [[70.0, 100.0]] * 5, 'RMS Noise': [[0.001, 0.002]] * 5}),
        'y_both_minerr': np.array([0.0001, 0.0001]), 'y_both_maxerr': np.array([0.0001, 0.0001]),
    }
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


@pytest.mark.parametrize('kind', ['histogram', 'scatter'])
@pytest.mark.parametrize('typeplotback', ['white', 'dark'])
@pytest.mark.parametrize('save', [False, True])
@pytest.mark.parametrize('inherit', [False, True])
def test_plot_functions_honor_requested_theme(plot_functions, monkeypatch, tmp_path, kind, typeplotback, save, inherit):
    background, foreground = plot_background_colors(typeplotback)
    displayed = []

    def inspect_display():
        figure = plt.gcf()
        displayed.append(figure)
        assert figure.get_facecolor() == to_rgba(background)
        for axis in figure.axes:
            assert axis.get_facecolor() == to_rgba(background)
            assert axis.xaxis.label.get_color() == foreground
            assert axis.yaxis.label.get_color() == foreground
            assert all(spine.get_edgecolor() == to_rgba(foreground) for spine in axis.spines.values())
        FigureCanvasAgg(figure).draw()
        np.testing.assert_array_equal(np.asarray(figure.canvas.buffer_rgba())[0, 0],
                                      np.array(to_rgba(background)) * 255)

    monkeypatch.setattr(plt, 'show', inspect_display)
    monkeypatch.chdir(tmp_path)
    plot_functions['plot_background_colors'] = plot_background_colors
    options = {} if inherit else {'typeplotback': typeplotback}
    outer_theme = typeplotback if inherit else ('dark' if typeplotback == 'white' else 'white')
    with plot_style(outer_theme):
        original = dict(plt.rcParams)
        try:
            if kind == 'histogram':
                assert plot_functions['outcome_dist_hists'](
                    [[0.1, 1.0, 10.0]] * 3, 'both', [1.0, 10.0, 100.0],
                    save=save, **options,
                ) is None
            else:
                data = pd.DataFrame({'Names': ['A', 'B'], 'l_mag': ['A', 'B'], 'l_kin': ['A', 'B'],
                                     'l_both': ['A', 'B'], 'x': [20.0, 40.0], 'y_both': [0.001, 0.01],
                                     'd': [20.0, 40.0], 's': [-2.0, -1.0]})
                before = data.copy(deep=True)
                errors = [[np.full(2, 0.0001), np.full(2, 0.0001)] for _ in range(3)]
                detections = pd.DataFrame({'Name': ['A', 'B'], 'Freq': [20.0, 40.0], 'Flux': [0.001, 0.01]})
                assert plot_functions['scatter_plot'](
                    data, 'both', errors, [np.ones(2), np.ones(2)], detections,
                    [1.0] * 4, zoom=True, save=save, **options,
                ) is None
                pd.testing.assert_frame_equal(data, before)
                np.testing.assert_array_equal(displayed[0].axes[0].collections[0].get_offsets(), data[['x', 'y_both']])
            assert dict(plt.rcParams) == original
            assert len(displayed) == 1
            if save:
                filename = tmp_path / ('hist.pdf' if kind == 'histogram' else 'zoom.pdf')
                assert filename.read_bytes().startswith(b'%PDF-')
        finally:
            plt.close('all')