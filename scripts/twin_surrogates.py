    
    
    # array object and fast numerics
import numpy as np
from numpy import random
import sys
#sys.path.append(".")
import pyximport
pyximport.install()#._ext/
#import direnv


#from numerics import _twins_s, _twin_surrogates
#from numerics import  
print("In module products sys.path[0], __package__ ==", sys.path[0], __package__)
#from ..utils import progressbar
from numerics import _twins_s, _twin_surrogates
#from ._ext import numerics
# easy progress bar handling


    
class Surrogates:

    """
    Encapsulates structures and methods related to surrogate time series.
    Provides data structures and methods to generate surrogate data sets from a
    set of time series and to evaluate the significance of various correlation
    measures using these surrogates.
    More information on time series surrogates can be found in [Schreiber2000]_
    and [Kantz2006]_.
    """    
    
    # FIXME: I(wb) included the line
    # dimension = embedding_array.shape[2]
    # whose missing caused an error. I can't guarantee if it is correct.

    def embed_time_series_array(self, time_series_array, dimension, delay):
            """
            Return a :index:`delay embedding` of all time series.
            .. note::
            Only works for scalar time series!
            **Example:**
            >>> ts = Surrogates.SmallTestData().original_data
            >>> Surrogates.SmallTestData().embed_time_series_array(
            ...     time_series_array=ts, dimension=3, delay=2)[0,:6,:]
            array([[ 0.        ,  0.61464833,  1.14988147],
                [ 0.31244015,  0.89680225,  1.3660254 ],
                [ 0.61464833,  1.14988147,  1.53884177],
                [ 0.89680225,  1.3660254 ,  1.6636525 ],
                [ 1.14988147,  1.53884177,  1.73766672],
                [ 1.3660254 ,  1.6636525 ,  1.76007351]])

            :type time_series_array: 2D array [index, time]
            :arg time_series_array: The time series array to be normalized.
            :arg int dimension: The embedding dimension.
            :arg int delay: The embedding delay.
            :rtype: 3D array [index, time, dimension]
            :return: the embedded time series.
            """
            if self.silence_level <= 1:
                print(f"Embedding all time series in dimension {dimension} "
                    f"and with lag {delay} ...")
            (N, n_time) = time_series_array.shape

            embedding = np.empty((N, n_time - (dimension - 1)*delay, dimension))


            return embedding


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