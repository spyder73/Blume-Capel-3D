#include <iostream>
#include <vector>
#include "random.hpp"
#include "triple.hpp"
#pragma once
using namespace std;

class BlumeCapel {
private:
    const double kT;
    const int L; //CUBE SIZE
    const double D;
    const double h;
    const double K=1.0/kT;
    triple<pair<const bool,const bool>,pair<const bool,const bool>,pair<const bool,const bool>> Boundaries;
    vector<vector<vector<int>>> lattice; //CUBE TENSOR
    int spin_i, spin_j, spin_k; //INDICES WHICH WILL BE UPDATED IN EACH ITERATION
    int* current_adress; //ADRESS OF THE SELECTED SPIN, UPDATED IN EACH ITERATION

    //Necessary for Calculation of Energy etc. depending on if periodic boundary conditions are ON/OFF:
    int boundary_i0, boundary_i1, boundary_j0, boundary_j1, boundary_k0, boundary_k1;
    int m_i0, m_i1, m_j0, m_j1, m_k0, m_k1;
public:

    //FUNCTION TO CHOOSE A NEW RANDOM SPIN ON THE CUBE LATTICE
    void randomizeSpins();

    //CONSTRUCTOR OF THIS GIVEN CLASS; CREATES RANDOMIZED CUBE TENSOR AND CHOOSES A RANDOM FIRST SPIN 
    BlumeCapel(const int L, const double kT, const double D, const double h, triple<pair<const bool,const bool>,pair<const bool,const bool>,pair<const bool,const bool>> Boundaries);

    ~BlumeCapel();

    double get_kT();

    //FUNCTION TO PRINT THE LATTICE
    void print();

    //RETURNS THE CURRENT ENERGY OF THE LATTICE
    double getE();

    //RETURNS THE POSSIBLE CHANGE IN ENERGY IF A SPIN IS TO BE FLIPPED
    double getdE();

    //RETURNS CHANGE IN ENERGY FOR SPIN-FLIP TO/FROM ZERO
    double getdE_zero(string b, string c);

    //TWO-POINT-CORRELATION FUNCTION
    double tp_correlation(triple<int,int,int>& spin_1, triple<int,int,int>& spin_2);
    
    //ALGORITHM FOR FLIPPING A RANDOM SPIN, BASED ON CALCULATION OF dE
    void spinFlip();

    //CALCULATES THE MAGNETIZATION AND SUSCEPTIBILITY AS A PAIR
    pair<double,double> calcMagChi();

    //WOLFF-ALGORITHM FOR CLUSTER FLIPS
    void clusterFlip();
};
