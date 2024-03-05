#include <iostream>
#include <vector>
#include "triple.hpp"
#pragma once
using namespace std;

//CALCULATE MEAN VALUE OF A GIVEN ARRAY WITH N VALUES
template<typename T>
double CalcMean(const vector<T>& k);

template<typename T>
double CalcMean_pow(const vector<T>& k, int power);

template<typename T>
triple<vector<T>, vector<T>, vector<T>> CalcBins_Mag(const vector<vector<T>>& matrix);

//CREATE BINS OF EQUALLY SIZED DATASETS
template<typename T>
vector<vector<T>> create_bins(const vector<T>& data, const int& N_bin, const int& bin_size);

//Calculate JackKnife Average Vector for any Distribution of Bins
template<typename T>
vector<T> CalcBins(const vector<vector<T>>& matrix);

template<typename T>
pair<T,T> U4_calc(const vector<T>& M2, const vector<T>& M4);

//CALCULATE JACKKNIFE ESTIMATE OF A GIVEN ARRAY, POWER USED FOR <X^POWER> AVERAGE, RETURNS VECTOR WITH N ESTIMATES FOR EACH i != j
template<typename T>
vector<T> jack_avg(const vector<T>& data, const int power);

//CALCULATE MEAN OF JACKKNIFE ESTIMATES + STD
template<typename T>
pair<T,T> jack_mean_std(const vector<T>& data);