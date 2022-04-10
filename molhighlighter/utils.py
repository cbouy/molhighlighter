from functools import wraps
from collections import namedtuple
from numpy import linspace
from colorsys import hls_to_rgb


Substitution = namedtuple("Substitution", ["content", "start", "end"])


def requires_config(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        if not self._is_configured:
            self.configure()
        return method(self, *args, **kwargs)
    return wrapper


sequential_palette = ["#93e467", "#e36262", "#62d4e3", "#6275e3", "#d262e3"]


def get_auto_palette(n_colors, h=.0, s=.75, l=.65):
    """Returns a perceptually uniform palette of colors
    
    Parameters
    ----------
    n_colors : int
        Number of colors to output
    h : float
        First color in the palette [0, 1]
    s : float
        Saturation [0, 1]
    l : float
        Lightness [0, 1]
    """
    hues = linspace(0, 1, n_colors + 1)[:-1]
    hues += h
    hues %= 1
    return [hls_to_rgb(h, l, s) for h in hues]
