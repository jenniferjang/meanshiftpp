import numpy as np
cimport numpy as np
from sklearn.neighbors import KDTree
from scipy.sparse import csr_matrix
from datetime import datetime
import pprint


cdef extern from "meanshiftpp.h":
    void generate_offsets_cy(int d, 
                             int base,
                             int * offsets)

    void shift_cy(int n,
                  int d,
                  int base,
                  float bandwidth,
                  int * offsets,
                  float * X_shifted)

cdef generate_offsets_np(d, 
                         base,
                         np.ndarray[np.int32_t, ndim=2, mode="c"] offsets):
    generate_offsets_cy(d, 
                        base,
                        <int *> np.PyArray_DATA(offsets))

cdef shift_np(n,
              d,
              base,
              bandwidth,
              np.ndarray[np.int32_t, ndim=2, mode="c"] offsets,
              np.ndarray[float, ndim=2, mode="c"] X_shifted):
    shift_cy(n,
             d,
             base,
             bandwidth,
             <int *> np.PyArray_DATA(offsets),
             <float *> np.PyArray_DATA(X_shifted))


class MeanshiftPP:
    """
    Parameters
    ----------
    
    bandwidth: Radius for binning points. Points are assigned to the bin 
               corresponding to floor division by bandwidth

    threshold: Stop shifting if the L2 norm between iterations is less than
               threshold

    iterations: Maximum number of iterations to run

    """

    def __init__(self, bandwidth, threshold=0.0001, iterations=None):
        self.bandwidth = bandwidth
        self.threshold = threshold
        self.iterations = iterations


    def fit_predict(self, X):
        """
        Determines the clusters in either `iterations` or when the L2 
        norm of consecutive iterations is less than `threshold`, whichever 
        comes first.
        Each shift has two steps: First, points are binned based on floor 
        division by bandwidth. Second, each bin is shifted to the 
        weighted mean of its 3**d neighbors. 
        Lastly, points that are in the same bin are clustered together.

        Parameters
        ----------
        X: Data matrix. Each row should represent a datapoint in 
           Euclidean space

        Returns
        ----------
        (n, ) cluster labels
        """
        
        X = np.ascontiguousarray(X, dtype=np.float32)
        n, d = X.shape
        X_shifted = np.copy(X)

        result = np.full(n, -1, dtype=np.int32)

        iteration = 0
        base =  3
        offsets = np.full((base**d, d), -1, dtype=np.int32)
        generate_offsets_np(d, base, offsets)
        
        while iteration < iterations or not iterations:
          #print("Iteration: ", cnt, len(np.unique(X, axis=0)))
          iteration += 1

          shift_np(n, d, base, self.bandwidth, offsets, X_shifted)

          if np.linalg.norm(np.subtract(X, X_shifted)) <= self.threshold:
            break

          X = np.copy(X_shifted)
          
        _, result = np.unique(X_shifted, return_inverse=True, axis=0)
        
        return result
