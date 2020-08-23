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
        time1 = datetime.now()
        generate_offsets_np(d, base, offsets)
        time2 = datetime.now()
        print("A", (time2-time1).total_seconds())
        
        while True:
          #print("Iteration: ", cnt, len(np.unique(X, axis=0)))
          cnt += 1

          time3 = datetime.now()
          shift_np(n, d, base, self.bandwidth, offsets, X_shifted)
          time4 = datetime.now()
          print("B", (time4-time3).total_seconds())

          if np.linalg.norm(np.subtract(X, X_shifted)) <= self.threshold:
            break

          time5 = datetime.now()
          print("C", (time5-time4).total_seconds())

          X = np.copy(X_shifted)
          time6 = datetime.now()
          print("D", (time6-time5).total_seconds())

        time7 = datetime.now()

        _, result = np.unique(X_shifted, return_inverse=True, axis=0)
        time8 = datetime.now()
        print("E", (time8-time7).total_seconds())
        
        return result
