import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.use('Agg')


def threshold_otsu(hist):
    """Return threshold value based on Otsu's method.

    hist : array, or 2-tuple of arrays, optional
        Histogram from which to determine the threshold, and optionally a
        corresponding array of bin center intensities.
        An alternative use of this function is to pass it only hist.
    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.
    References
    ----------
    .. [1] Wikipedia, https://en.wikipedia.org/wiki/Otsu's_Method
    """
    counts, bin_centers = hist
    bin_centers = bin_centers[:-1]

    # class probabilities for all possible thresholds
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(counts * bin_centers) / weight1
    mean2 = (np.cumsum((counts * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of ``weight1``/``mean1`` should pair with zero values in
    # ``weight2``/``mean2``, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    if len(variance12) == 0:
        return 0
    idx = np.nanargmax(variance12)
    threshold = bin_centers[idx]

    return threshold


def array2hist(array, binWidth=0.2):
    counts, bins = np.histogram(array, bins=np.arange(0, max(array)+binWidth, binWidth))
    return counts, bins


def makePlot(hist, thresh, fname):
    counts, bins = hist
    plt.hist(bins[:-1], bins, weights=counts)
    plt.axvline(thresh, color='r')
    plt.savefig(fname)
    plt.close()
