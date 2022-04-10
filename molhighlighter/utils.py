import numpy as np
from colorsys import hls_to_rgb

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
    hues = np.linspace(0, 1, n_colors + 1)[:-1]
    hues += h
    hues %= 1
    return [hls_to_rgb(h, l, s) for h in hues]
