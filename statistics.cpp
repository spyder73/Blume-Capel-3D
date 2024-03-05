#include <iostream>
#include <vector>
#include "statistics.hpp"
using namespace std;

//CALCULATE MEAN VALUE OF A GIVEN ARRAY WITH N VALUES
template<typename T>
double CalcMean(const vector<T>& k) {
    double sum = 0.0;
    for (auto i = k.begin(); i != k.end(); i++) {
        sum += *i;
    }
    return sum/static_cast<double>(k.size());
}

template<typename T>
double CalcMean_pow(const vector<T>& k, int power) {
    double sum = 0.0;
    for (auto i = k.begin(); i != k.end(); i++) {
        sum += pow(*i,power);
    }
    return sum/static_cast<double>(k.size());
}

//CREATE BINS OF EQUALLY SIZED DATASETS
template<typename T>
vector<vector<T>> create_bins(const vector<T>& data, const int& N_bin, const int& bin_size) {
    vector<vector<T>> matrix(N_bin, vector<T>(bin_size));
    size_t k=0;
    for(int i=0; i<N_bin; i++) {
        for (int j=0; j<bin_size; j++) {
            matrix[i][j] = data[k];
            k++;
        }
    }
    return matrix;
}

//Calculate All Binned Data (MEAN of each Bin), POWER = 1 (FIXED)
template<typename T>
vector<T> CalcBins(const vector<vector<T>>& matrix) {
    vector<T> M;
    for(auto i = matrix.begin(); i != matrix.end(); i++) {
        vector<T> temp;
        for(auto j = i->begin(); j != i->end(); j++) {
            temp.push_back(*j);
        }
        M.push_back(CalcMean(temp));
    }
    return M;
}

//Calculate Binned Data with power for Magnetization:
template<typename T>
triple<vector<T>, vector<T>, vector<T>> CalcBins_Mag(const vector<vector<T>>& matrix) {
    vector<T> M1;
    vector<T> M2;
    vector<T> M4;
    for(auto i = matrix.begin(); i != matrix.end(); i++) {
        vector<T> temp;
        for(auto j = i->begin(); j != i->end(); j++) {
            temp.push_back(*j);
        }
        M1.push_back(CalcMean_pow(temp, 1));
        M2.push_back(CalcMean_pow(temp, 2));
        M4.push_back(CalcMean_pow(temp, 4));
    }
    return triple<vector<T>, vector<T>, vector<T>>(M1,M2,M4);
}


template<typename T>
pair<T,T> U4_calc(const vector<T>& M2, const vector<T>& M4) {
    vector<T> f_m;
    for(size_t i=0; i<M4.size(); i++) {
        f_m.push_back(M4[i]/pow(M2[i],2));
    }
    pair<T,T> result = jack_mean_std(f_m);
    return pair<T,T>(result.first, result.second);
}

//CALCULATE JACKKNIFE ESTIMATE OF A GIVEN ARRAY, POWER USED FOR <X^POWER> AVERAGE, RETURNS VECTOR WITH N ESTIMATES FOR EACH i != j
template<typename T>
vector<T> jack_avg(const vector<T>& data, const int power) {
    vector<T> y;
    double y_i = 1.0/(data.size()-1.0);
    for (auto i=data.begin(); i!=data.end(); i++) {
        T sum = 0.0;
        for (auto j=data.begin(); j!=data.end(); j++) {
            if(j!=i) {
                sum += pow(*j,power);
            }
        }
        y.push_back(y_i*sum);
    }
    return y;
}

//CALCULATE MEAN OF JACKKNIFE ESTIMATES + STD
template<typename T>
pair<T,T> jack_mean_std(const vector<T>& data) {
    T mean = CalcMean(data);
    T val = 0.0;
    for (auto i = data.begin(); i != data.end(); i++) {
        val += pow(*i - mean, 2);
    }
    T jack_std = sqrt(data.size()-1)/sqrt(data.size()) * sqrt(val);
    return pair<T,T>(mean, jack_std);
}


// Explicit instantiations
template double CalcMean(const vector<double>& k);
template double CalcMean_pow(const vector<double>&k, int power);
template vector<vector<double>> create_bins(const vector<double>& data, const int& N_bin, const int& bin_size);
template vector<double> CalcBins(const vector<vector<double>>& matrix);
template pair<double,double> U4_calc(const vector<double>& M2, const vector<double>& M4);
template vector<double> jack_avg(const vector<double>& data, const int power);
template pair<double,double> jack_mean_std(const vector<double>& data);

template triple<vector<double>, vector<double>, vector<double>> CalcBins_Mag(const vector<vector<double>>& matrix);