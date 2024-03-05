#include <iostream>
#include <vector>
#include "blume.hpp"
using namespace std;


//FUNCTION TO CHOOSE A NEW RANDOM SPIN ON THE CUBE LATTICE
void BlumeCapel::randomizeSpins() {
    spin_i = rdm_int(0, static_cast<int>(L)-1);
    spin_j = rdm_int(0, static_cast<int>(L)-1);
    spin_k = rdm_int(0, static_cast<int>(L)-1);
    current_adress = &lattice[spin_i][spin_j][spin_k];

    if (spin_i == 0 && !Boundaries.first.first) {boundary_i0 = 0; m_i0 = 0;}
    else {boundary_i0 = (spin_i - 1 + L )%L; m_i0 = 1;}

    if (spin_i == (L-1) && !Boundaries.first.second) {boundary_i1 = 0; m_i1 = 0;}
    else {boundary_i1 = (spin_i + 1)%L; m_i1 = 1; }

    if (spin_j == 0 && !Boundaries.second.first) {boundary_j0 = 0; m_j0 = 0;}
    else {boundary_j0 = (spin_j - 1 + L )%L; m_j0 = 1;}

    if (spin_j == (L-1) && !Boundaries.second.second) {boundary_j1 = 0; m_j1 = 0;}
    else {boundary_j1 = (spin_j + 1)%L; m_j1 = 1;}
     
    if (spin_k == 0 && !Boundaries.third.first) {boundary_k0 = 0; m_k0 = 0;}
    else {boundary_k0 = (spin_k - 1 + L )%L; m_k0 = 1;}

    if (spin_k == (L-1) && !Boundaries.third.second) {boundary_k1 = 0; m_k1 = 0;}
    else {boundary_k1 = (spin_k + 1)%L; m_k1 = 1;}
}

//CONSTRUCTOR OF THIS GIVEN CLASS; CREATES RANDOMIZED CUBE TENSOR AND CHOOSES A RANDOM FIRST SPIN
//Update Blume-Capel-Parameters to be dimensionless (D/ß & h/ß):
BlumeCapel::BlumeCapel(const int L, const double kT, const double D, const double h, triple<pair<const bool,const bool>,pair<const bool,const bool>,pair<const bool,const bool>> Boundaries): 
    kT(kT), L(L), D(D*kT), h(h*kT), Boundaries(Boundaries), lattice(L, vector<vector<int>>(L, vector<int>(L,0))) {

    for(auto i=lattice.begin(); i!=lattice.end(); i++) {
        for(auto j=i->begin(); j!=i->end(); j++) {
            for(auto k=j->begin(); k!=j->end(); k++) {
                *k = rdm_int(-1,1);
            }
        }
    }
    randomizeSpins();

}

BlumeCapel::~BlumeCapel() {}

double BlumeCapel::get_kT() {
    return kT;
}

//FUNCTION TO PRINT THE LATTICE
void BlumeCapel::print() {
for(auto i=lattice.begin(); i!=lattice.end(); i++) {
    for(auto j=i->begin(); j!=i->end(); j++) {
        for(auto k=j->begin(); k!=j->end(); k++) {
            cout << *k << " ";
            }
    cout << endl;
        }
    cout << endl << endl;
    }
}

//RETURNS THE CURRENT ENERGY OF THE LATTICE
double BlumeCapel::getE() {
    double E = 0.0;
    for (int i = 0; i != L; i++) {
        for (int j = 0; j != L; j++) {
            for (int k=0; k != L; k++) {
                if(i != 0 || (i == 0 && Boundaries.first.first)) {E += lattice[i][j][k] * (lattice[(i-1+L) % L][j][k]);}
                if(i != (L-1) || (i==(L-1) && Boundaries.first.second)) {E += lattice[i][j][k] * (lattice[(i+1) % L][j][k]);}
                if(j != 0 || (j == 0 && Boundaries.second.first)) {E += lattice[i][j][k] * (lattice[i][(j-1+L) % L][k]);}
                if(j != (L-1) || (j==(L-1) && Boundaries.second.second)) {E += lattice[i][j][k] * (lattice[i][(j+1) % L][k]);}
                if(k != 0 || (k == 0 && Boundaries.third.first)) {E += lattice[i][j][k] * (lattice[i][j][(k-1+L) % L]);}
                if(k != (L-1) || (k == (L-1) && Boundaries.third.second)) {E += lattice[i][j][k] * (lattice[i][j][(k+1) % L]);}
                E += D*pow(lattice[i][j][k],2) - h*lattice[i][j][k];
                }
            }
        }
    return E;
}

//RETURNS THE POSSIBLE CHANGE IN ENERGY IF A SPIN IS TO BE FLIPPED
double BlumeCapel::getdE() {
    double dE=0.0;
    dE += 2*lattice[spin_i][spin_j][spin_k] * (lattice[boundary_i0][spin_j][spin_k]*m_i0 + lattice[boundary_i1][spin_j][spin_k]*m_i1);
    dE += 2*lattice[spin_i][spin_j][spin_k] * (lattice[spin_i][boundary_j0][spin_k]*m_j0 + lattice[spin_i][boundary_j1][spin_k]*m_j1);
    dE += 2*lattice[spin_i][spin_j][spin_k] * (lattice[spin_i][spin_j][boundary_k0]*m_k0 + lattice[spin_i][spin_j][boundary_k1]*m_k1);
    dE += 2*h*lattice[spin_i][spin_j][spin_k];
    return dE;
}

//RETURNS CHANGE IN ENERGY FOR SPIN-FLIP TO/FROM ZERO
double BlumeCapel::getdE_zero(string b, string c) {
    double dE;
    if (b == "to") {
        double E_i = 0.0; double E_f = 0.0;
        E_i -= lattice[spin_i][spin_j][spin_k] * (lattice[boundary_i0][spin_j][spin_k]*m_i0 + lattice[boundary_i1][spin_j][spin_k]*m_i1);
        E_i -= lattice[spin_i][spin_j][spin_k] * (lattice[spin_i][boundary_j0][spin_k]*m_j0 + lattice[spin_i][boundary_j1][spin_k]*m_j1);
        E_i -= lattice[spin_i][spin_j][spin_k] * (lattice[spin_i][spin_j][boundary_k0]*m_k0 + lattice[spin_i][spin_j][boundary_k1]*m_k1);
        E_i += D - h*lattice[spin_i][spin_j][spin_k];
        dE = E_f - E_i;
    }
    if (b == "from") {
        double E_i = 0.0; double E_f = 0.0;
        if (c == "+") {
            E_f -= 1 * (lattice[boundary_i0][spin_j][spin_k]*m_i0 + lattice[boundary_i1][spin_j][spin_k]*m_i1);
            E_f -= 1 * (lattice[spin_i][boundary_j0][spin_k]*m_j0 + lattice[spin_i][boundary_j1][spin_k]*m_j1);
            E_f -= 1 * (lattice[spin_i][spin_j][boundary_k0]*m_k0 + lattice[spin_i][spin_j][boundary_k1]*m_k1);
            E_f += D - h; // + D*(1)^2 - h * (1)
        }
        if (c == "-") {
            E_f += 1 * (lattice[boundary_i0][spin_j][spin_k]*m_i0 + lattice[boundary_i1][spin_j][spin_k]*m_i1);
            E_f += 1 * (lattice[spin_i][boundary_j0][spin_k]*m_j0 + lattice[spin_i][boundary_j1][spin_k]*m_j1);
            E_f += 1 * (lattice[spin_i][spin_j][boundary_k0]*m_k0 + lattice[spin_i][spin_j][boundary_k1]*m_k1);
            E_f += D + h; // + D*(-1)^2 - h *(-1)
        }
        dE = E_f - E_i;
    }
    return dE;
}

//ALGORITHM FOR FLIPPING A RANDOM SPIN, BASED ON CALCULATION OF dE
void BlumeCapel::spinFlip() {
    //SPIN = 0 CASE:
    if(lattice[spin_i][spin_j][spin_k] == 0) {
        double r1 = rdm_double(0.0, 1.0);
        if (r1 <= 0.5) {
            double dE = getdE_zero("from", "+");
            double f = rdm_double(0.0, 1.0);
            if (dE >= 0) {
                if (f < exp(-dE/kT)) {
                    *current_adress = 1;
                }
            }
            else {
                *current_adress = 1;
            }
        randomizeSpins();
        }

        else if (r1 > 0.5) {
            double dE = getdE_zero("from", "-");
            double f = rdm_double(0.0, 1.0);
            if (dE >= 0) {
                if (f < exp(-dE/kT)) {
                    *current_adress = -1;
                }
            }
            else {
                *current_adress = -1;
            }
        randomizeSpins();
        }
    }
    if(lattice[spin_i][spin_j][spin_k] != 0) {
        double r2 = rdm_double(0.0, 1.0);
        if (r2 <= 0.5) {
            double dE = getdE_zero("to", "0");
            double f = rdm_double(0.0, 1.0);
            if (dE >= 0) {
                if (f < exp(-dE/kT)) {
                    *current_adress = 0;
                }
            }
            else {
                *current_adress = 0;
            }
            randomizeSpins();
        }
        else if (r2 > 0.5) {
            double dE = getdE();
            double f = rdm_double(0.0, 1.0);
            if (dE >= 0) {
                if (f < exp(-dE/kT)) {
                    *current_adress *= -1;
                }
            }
            else {
                *current_adress *= -1;
            }
            randomizeSpins();
        }
    }
}

//CALCULATES THE MAGNETIZATION AND SUSCEPTIBILITY AS A PAIR
pair<double,double> BlumeCapel::calcMagChi() {
    double M=0.0;
    for (auto i=lattice.begin(); i!= lattice.end(); i++) {
        for (auto j=i->begin(); j!= i->end(); j++) {
            for(auto k=j->begin(); k != j->end(); k++) {
                M += *k;
            }
        }
    }
    return pair<double,double>(static_cast<double>(abs(M)), static_cast<double>(pow(M,2))/(L*L*L));
}

//TWO-POINT-CORRELATION FUNCTION
double BlumeCapel::tp_correlation(triple<int,int,int>& spin_1, triple<int,int,int>& spin_2) {
    int i = spin_1.first; int j = spin_1.second; int k = spin_1.third;
    int l = spin_2.first; int m = spin_2.second; int n = spin_2.third;
    return lattice[(i+L)%L][(j+L)%L][(k+L)%L] * lattice[(l+L)%L][(m+L)%L][(n+L)%L];
}

//WOLFF-ALGORITHM FOR CLUSTER FLIPS
void BlumeCapel::clusterFlip() {
    vector<vector<vector<bool>>> visited(L, vector<vector<bool>>(L, vector<bool>(L, false))); //No Spin is visited in the beginning

    int target = lattice[spin_i][spin_j][spin_k];

    if (target != 0) {
        //Create Stack
        vector< triple<int,int,int> > stack;
        stack.push_back({spin_i,spin_j,spin_k});

        while (!stack.empty()) {
            int i = stack.back().first;
            int j = stack.back().second;
            int k = stack.back().third;

            stack.pop_back();

            // Check if visited already
            if (!visited[i][j][k]) {
                visited[i][j][k] = true;
                lattice[i][j][k] *= -1;

            if( (i==0 && Boundaries.first.first) || i!=0 ) {
                    if (lattice[(i - 1 + L) % L][j][k] == target) {
                        // Add neighboring spins to the stack
                        double f = rdm_double(0.0, 1.0);
                        if (f < (1.0 - exp(-2.0 * K))) {
                            stack.push_back({(i - 1 + L) % L, j, k});
                        }
                    }
                }

                if ( (i == (L-1) && Boundaries.first.second) || i != (L-1) ) {
                    if (lattice[(i + 1) % L][j][k] == target) {
                        // Add neighboring spins to the stack
                        double f = rdm_double(0.0, 1.0);
                        if (f < (1.0 - exp(-2.0 * K))) {
                            stack.push_back({(i + 1) % L, j, k});
                        }
                    }
                }

                if ( (j == 0 && Boundaries.second.first) || j != 0 ) {
                    if (lattice[i][(j - 1 + L) % L][k] == target) {
                        // Add neighboring spins to the stack
                        double f = rdm_double(0.0, 1.0);
                        if (f < (1.0 - exp(-2.0 * K))) {
                            stack.push_back({i, (j - 1 + L) % L, k});
                        }
                    }
                }

                if ( (j == (L-1) && Boundaries.second.second) || j != (L-1) ) {
                    if (lattice[i][(j + 1) % L][k] == target) {
                        // Add neighboring spins to the stack
                        double f = rdm_double(0.0, 1.0);
                        if (f < (1.0 - exp(-2.0 * K))) {
                            stack.push_back({i, (j + 1) % L, k});
                        }
                    }
                }

                if ( (k == 0 && Boundaries.third.first) || k != 0 ) {
                    if (lattice[i][j][(k - 1 + L) % L] == target) {
                        // Add neighboring spins to the stack
                        double f = rdm_double(0.0, 1.0);
                        if (f < (1.0 - exp(-2.0 * K))) {
                            stack.push_back({i, j, (k - 1 + L) % L});
                        }
                    }
                }

                if ( (k == (L-1) && Boundaries.third.second) || k != (L-1) ) {
                    if (lattice[i][j][(k + 1) % L] == target) {
                        // Add neighboring spins to the stack
                        double f = rdm_double(0.0, 1.0);
                        if (f < (1.0 - exp(-2.0 * K))) {
                            stack.push_back({i, j, (k + 1) % L});
                        }
                    }
                }
            }
        }
        randomizeSpins();
    }
    else {randomizeSpins();}
}