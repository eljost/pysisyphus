from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


def fig_ax_from_canvas_agg():
    """Figure and axis w/o tkinter issue: RuntimeError main thread is not in main loop."""
    fig = Figure()
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot()
    return fig, ax
