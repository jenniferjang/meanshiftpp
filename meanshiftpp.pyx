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
    
    bandwidth: Number of neighbors required for a point to be labeled a core point. Works
               in conjunction with eps_density

    """

    def __init__(self, bandwidth, threshold=0.01):
        self.bandwidth = bandwidth
        self.threshold = threshold


    def fit_predict(self, X):
        """
        Parameters
        ----------
        

        Returns
        ----------
        (n, ) cluster labels
        """
        
        X = np.ascontiguousarray(X, dtype=np.float32)
        n, d = X.shape
        X_shifted = np.copy(X)

        result = np.full(n, -1, dtype=np.int32)

        cnt = 0
        base =  3
        offsets = np.full((base**d, d), -1, dtype=np.int32)
        generate_offsets_np(d, base, offsets)
        
        while True:
          #print("Iteration: ", cnt, len(np.unique(X, axis=0)))
          cnt += 1

          shift_np(n, d, base, self.bandwidth, offsets, X_shifted)
          
          if np.linalg.norm(np.subtract(X, X_shifted)) <= self.threshold:
            break

          X = np.copy(X_shifted)
          
        _, result = np.unique(X_shifted, return_inverse=True, axis=0)
        
        return result
