from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


def fig_ax_from_canvas_agg():
    fig = Figure()
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot()
    return fig, ax
