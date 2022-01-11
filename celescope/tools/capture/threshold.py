import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.use('Agg')

# minimum array length to use otsu method
OTSU_MIN_LEN = 50

class Otsu():
    """
    Return 1 if len(array) < OTSU_MIN_LEN
    """
    def __init__(self, array, log_base=10):
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


    def _threshold_otsu(self):
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

    def _array2hist(self, binWidth=0.2):
        self.counts, self.bins = np.histogram(self.array, bins=np.arange(0, max(self.array)+binWidth, binWidth))

    def make_plot(self, plot_path):
        plt.hist(x=self.bins[:-1], bins=self.bins, weights=self.counts)
        plt.axvline(self.threshold, color='r')
        plt.xlabel(f'log{self.log_base} observed read/UMI counts')
        plt.ylabel('Frequency')
        plt.savefig(plot_path)
        plt.close()

    def run(self):
        """
        return threshold
        """
        if len(self.array) < OTSU_MIN_LEN:
            return 1
        self._array2hist()
        self._threshold_otsu()
        return_threshold = int(self.log_base ** self.threshold)

        return return_threshold


class Auto():
    """
    threshold = top 1% cell count / 10
    """
    def __init__(self, array):
        self.array = array
    
    def run(self):
        array = self.array
        if not array:
            return 1
        n_cell_1_percentile = len(array) // 100
        sorted_counts = sorted(array, reverse=True)
        count_cell_1_percentile = sorted_counts[n_cell_1_percentile]
        threshold = int(count_cell_1_percentile / 10)

        return threshold


class Threshold():
    """
    Args:
        array: array-like
        threshold_method: ['otsu', 'auto', 'hard']
        otsu_plot_path: str
        hard_threshold: int
    """
    def __init__(self, array, threshold_method='auto', otsu_plot_path=None, hard_threshold=None):
        self.array = [x for x in array if x > 0 ]
        self.threshold_method = threshold_method
        self.otsu_plot_path = otsu_plot_path
        self.hard_threshold = hard_threshold
    
    def run(self):
        """
        return threshold
        """
        if not self.array:
            return 1
        if self.threshold_method == 'otsu':
            otsu = Otsu(self.array)
            threshold = otsu.run()
            if self.otsu_plot_path:
                otsu.make_plot(self.otsu_plot_path)
        elif self.threshold_method == 'auto':
            auto = Auto(self.array)
            threshold = auto.run()
        elif self.threshold_method == 'hard':
            if self.hard_threshold:
                threshold = int(self.hard_threshold)
            else:
                raise Exception('hard_threshold must be set')
        else:
            raise ValueError(f'Unknown threshold method: {self.threshold_method}')

        return threshold