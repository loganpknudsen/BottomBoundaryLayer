"""Colour-scale limits computed from the data being plotted."""
import numpy as np


def _nice_up(x):
    """Smallest 1, 2 or 5 x 10^n that is >= x."""
    e = np.floor(np.log10(x))
    for m in (1, 2, 5, 10):
        if m * 10**e >= x:
            return float(m * 10**e)


def symlog_limits(panels):
    """Symmetric-log colour limits for a set of panels (array-likes).

    vmax:       the largest |value| over all panels, rounded up to 1, 2 or 5 x 10^n;
    lin_thresh: the smallest per-panel max |value| (the weakest panel), rounded down to a power
                of ten, so that panel just reaches the edge of the linear region;
    ticks:      0 and +/- each power of ten from lin_thresh up to vmax.
    """
    maxima = [float(np.nanmax(np.abs(np.asarray(p)))) for p in panels]
    vmax = _nice_up(max(maxima))
    lin_thresh = float(10 ** np.floor(np.log10(min(maxima))))
    decades = 10.0 ** np.arange(np.log10(lin_thresh), np.floor(np.log10(vmax)) + 1)
    ticks = [0.0] + [s * d for d in decades for s in (1, -1)]
    return vmax, lin_thresh, ticks
