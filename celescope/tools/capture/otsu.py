import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.use('Agg')

class Otsu():
    def __init__(self, array, plot_path, log_base=10):
        self.log_base = log_base
        if log_base == 2:
            self.array = np.log2(array)
        elif log_base == 10:
            self.array = np.log10(array)
        else:
            raise Exception('log_base must be 2, 10')

        self.threshold = 1
        self.counts = None
        self.bins = None

        # out
        self.plot_path = plot_path


    def threshold_otsu(self):
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
        counts = self.counts
        bin_centers = self.bins[:-1]

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
        self.threshold = bin_centers[idx]

    def array2hist(self, binWidth=0.2):
        self.counts, self.bins = np.histogram(self.array, bins=np.arange(0, max(self.array)+binWidth, binWidth))

    def make_plot(self):
        plt.hist(x=self.bins[:-1], bins=self.bins, weights=self.counts)
        plt.axvline(self.threshold, color='r')
        plt.xlabel(f'log{self.log_base} observed read/UMI counts')
        plt.ylabel('Frequency')
        plt.savefig(self.plot_path)
        plt.close()

    def run(self):
        """
        return threshold
        """
        self.array2hist()
        self.threshold_otsu()
        self.make_plot()
        return_threshold = int(self.log_base ** self.threshold)

        return return_threshold