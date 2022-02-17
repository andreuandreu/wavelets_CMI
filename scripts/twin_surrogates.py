    
    
    # array object and fast numerics
import numpy as np
from numpy import random

from ._ext.numerics import _embed_time_series_array, _recurrence_plot, \
    _twins_s, _twin_surrogates, _test_pearson_correlation, \
    _test_mutual_information

# easy progress bar handling
from ..utils import progressbar

    
    
    
# FIXME: I(wb) included the line
# dimension = embedding_array.shape[2]
# whose missing caused an error. I can't guarantee if it is correct.
def twins(self, embedding_array, threshold, min_dist=7):
    """
    Return list of the :index:`twins <pair: twins; surrogates>` of each
    state vector for all time series.
    Two state vectors are said to be twins if they share the same
    recurrences, i.e., if the corresponding rows or columns in the
    recurrence plot are identical.
    References: [Thiel2006]_, [Marwan2007]_.
    :type embedding_array: 3D array [index, time, dimension]
    :arg  embedding_array: The embedded time series array.
    :arg float threshold: The recurrence threshold.
    :arg number min_dist: The minimum temporal distance for twins.
    :rtype: [[number]]
    :return: the list of twins for each state vector in the time series.
    """
    if self.silence_level <= 1:
        print("Finding twins...")

    N = embedding_array.shape[0]
    n_time = embedding_array.shape[1]
    dimension = embedding_array.shape[2]
    twins = []

    #  Initialize the R matrix with ones
    R = np.empty((n_time, n_time))
    #  Initialize array to store the number of neighbors for each sample
    nR = np.empty(n_time)

    _twins_s(N, n_time, dimension, threshold, min_dist,
                embedding_array.copy(order='c'), R.copy(order='c'),
                nR.copy(order='c'), twins)

    return twins

#
#  Define methods to generate sets of surrogate time series
#


def twin_surrogates(self, original_data, dimension, delay, threshold,
                    min_dist=7):
    """
    Return surrogates using the twin surrogate method.
    Scalar twin surrogates are created by isolating the first component
    (dimension) of the twin surrogate trajectories.
    Twin surrogates share linear and nonlinear properties with the original
    time series, since they correspond to realizations of trajectories of
    the same dynamical systems with different initial conditions.
    References: [Thiel2006]_ [*], [Marwan2007]_.
    The twin lists of all time series are cached to facilitate a faster
    generation of several surrogates for each time series. Hence,
    :meth:`clear_cache` has to be called before generating twin surrogates
    from a different set of time series!
    :type original_data: 2D array [index, time]
    :arg original_data: The original time series.
    :arg int dimension: The embedding dimension.
    :arg int delay: The embedding delay.
    :arg float threshold: The recurrence threshold.
    :arg number min_dist: The minimum temporal distance for twins.
    :rtype: 2D array [index, time]
    :return: the twin surrogates.
    """
    #  The algorithm proceeds in several steps:
    #  1. Embed the original_data time series, using time delay embedding
    #     for simplicity. Use the same dimension and time delay delay for
    #     all time series for simplicity. Determine delay using time
    #     delayed mutual information and d using false nearest neighbors
    #     methods.
    #  2. Use the algorithm proposed in [*] to find twins
    #  3. Reconstruct one-dimensional twin surrogate time series

    (N, n_time) = original_data.shape
    n_time = n_time - (dimension-1)*delay

    #  Make sure that twins are calculated only once
    if self._twins_cached:
        twins = self._twins
    else:
        embedding = self.embed_time_series_array(original_data,
                                                    dimension, delay)
        twins = self.twins(embedding, threshold, min_dist)
        self._twins = twins
        self._twins_cached = True

    return _twin_surrogates(N, n_time, twins,
                            original_data.copy(order='c'))

#
#  Defines methods to generate correlation measure matrices based on
#  original_data and surrogate data for significance testing.
#