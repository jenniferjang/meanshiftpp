#include <iostream>
#include <unordered_map> 
#include <vector>
#include <string>
#include <sstream>
#include <chrono>

using namespace std;


void generate_offsets_cy(int d,
                         int base,
                         int * offsets) {
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


string _to_str(int d, 
               int * bin) {
    ostringstream os;
    for (int k = 0; k < d; k++) {
        os << bin[k] << ",";
    }
    return string(os.str());
}


void shift_cy(int n,
              int d,
              int base,
              float bandwidth,
              int * offsets, // 3**d by d
              float * X_shifted) {  // n by d
    /*
        
    */

    string bin_str;
    unordered_map< string, pair< vector<float>, int> > means;

    int * bins = new int[n * d];
    int * new_bin = new int[d];

    for (int i = 0; i < n; i++) {

        for (int k = 0; k < d; k++) {
            bins[i * d + k] = X_shifted[i * d + k] / bandwidth;
        }

        bin_str = _to_str(d, &bins[i * d]);

        if (means.find(bin_str) == means.end()) {
            means[bin_str] = make_pair(std::vector<float>(d, 0), 0);
        }
    }
    std::chrono::system_clock::time_point time1, time2, time3, time4, time5, time6;
    std::chrono::duration<double> elapsed_seconds;

    for (int i = 0; i < n; i++) {

        time1 = std::chrono::system_clock::now();
        for (int j = 0; j < pow(base, d); j++) {
            time2 = std::chrono::system_clock::now();
            for (int k = 0; k < d; k++) {
                 new_bin[k] = bins[i * d + k] + offsets[j * d + k];
            }

            time3 = std::chrono::system_clock::now();
            elapsed_seconds = time3-time2;
            cout << "a " << elapsed_seconds.count() << endl;

            bin_str = _to_str(d, new_bin);
            
            time4 = std::chrono::system_clock::now();
            elapsed_seconds = time4-time3;
            cout << "b " << elapsed_seconds.count() << endl;
            
            if (means.find(bin_str) != means.end()) {
                
                for (int k = 0; k < d; k++) {
                    means[bin_str].first[k] += X_shifted[i * d + k];
                }

                means[bin_str].second++;
            }
            
            time5 = std::chrono::system_clock::now();
            elapsed_seconds = time5-time4;
            cout << "c " << elapsed_seconds.count() << endl;
        }
        time6 = std::chrono::system_clock::now();
        elapsed_seconds = time6 - time1;
        cout << "d " << elapsed_seconds.count() << endl;
    }

    for (int i = 0; i < n; i++) {
        bin_str = _to_str(d, &bins[i * d]);
        
        for (int k = 0; k < d; k++) {
            X_shifted[i * d + k] = means[bin_str].first[k] * 1.0 / means[bin_str].second;
        }
    }
    
    delete[] bins;
    delete[] new_bin;
}
