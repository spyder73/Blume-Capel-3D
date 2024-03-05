#include <iostream>
#include <random>
#include <limits>
#include <utility>
#include "random.hpp"
using namespace std;

mt19937_64 generator;

int rdm_int(int a, int b) {
    uniform_int_distribution<int> distribution(a,b);
    int random_number = distribution(generator);
    return random_number;
}

double rdm_double(double a, double b) {
    uniform_real_distribution<double> distribution(a,b);
    double random_number = distribution(generator);
    return random_number;
}