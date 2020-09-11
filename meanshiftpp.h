#include <map> 
#include <vector>

using namespace std;


void generate_offsets_cy(int d,
                         int base,
                         int * offsets) {
    /*
        Generate 3**d neighbors for any point.

        Parameters
        ----------
        d: Dimensions
        base: 3, corresponding to (-1, 0, 1)
        offsets: (3**d, d) array of offsets to be added to 
                 a bin to get neighbors

    */

    int tmp_i;

    for (int i = 0; i < pow(base, d); i++) {
        tmp_i = i;
        for (int j = 0; j < d; j++) {
            if (tmp_i == 0) break;
            offsets[i * d + j] = tmp_i % base - 1;
            tmp_i /= base;
        }
    }
}


void shift_cy(int n,
              int d,
              int base,
              float bandwidth,
              int * offsets,
              float * X_shifted) {
    /*
        Generate 3**d neighbors for any point.

        Parameters
        ----------
        n: Number of samples
        d: Dimension
        base: 3, corresponding to (-1, 0, 1)
        bandwith: Radius for binning points. Points are assigned to the 
                  (d, ) bin corresponding to floor division by bandwidth
        offsets: (3**d, d) array of offsets to be added to 
                 a bin to get its neighbors
        X_shifted: (n, d) array of new points after one iteration
                   of shift

    */

    map< vector<int>, pair< vector<float>, int> > means;

    int * bins = new int[n * d];
    vector<int> current_bin(d);

    for (int i = 0; i < n; i++) {

        // Bin point
        for (int k = 0; k < d; k++) {
            bins[i * d + k] = X_shifted[i * d + k] / bandwidth;
            current_bin[k] = bins[i * d + k];
        }

        if (means.find(current_bin) == means.end()) {
            means[current_bin] = make_pair(std::vector<float>(d, 0), 0);
        }
    }
    
    for (int i = 0; i < n; i++) {

        for (int j = 0; j < pow(base, d); j++) {

            // Get neighbor
            for (int k = 0; k < d; k++) {
                 current_bin[k] = bins[i * d + k] + offsets[j * d + k];
            }

            // If neighbor exists, add it to the mean
            if (means.find(current_bin) != means.end()) {
                
                for (int k = 0; k < d; k++) {
                    means[current_bin].first[k] += X_shifted[i * d + k];
                }

                means[current_bin].second++;
            }
        }
    }

    // Set every point to the mean of its neighbors 
    for (int i = 0; i < n; i++) {

        for (int k = 0; k < d; k++) {
            current_bin[k] = bins[i * d + k];
        }
        
        for (int k = 0; k < d; k++) {
            X_shifted[i * d + k] = means[current_bin].first[k] * 1.0 / means[current_bin].second;
        }
    }
    
    delete[] bins;
}
