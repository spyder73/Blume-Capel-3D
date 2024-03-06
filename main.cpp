#include <cmath>
#include <fstream>
#include <chrono>
#include "blume.hpp"
#include "statistics.hpp"
using namespace std;

// Function definition
bool checkCinFail() {
    if (cin.fail()) {
        cerr << "Error: Input is invalid." << endl;
        cin.clear(); // Clear the error state of cin
        cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Ignore the rest of the input line
        return true; // Indicate that an error occurred
    }
    return false; // Indicate that input was successful
}

//High-Temperature-Susceptibility Expansion
double expansion(double temp) {
    double beta = 1.0/temp;
    double x = 0.0;

    x += 0.5308447611308816674 + 1.6907769625206168952*beta + 4.8619910172171602715*pow(beta,2)
    + 13.927143445449825611 * pow(beta, 3) + 38.903779467013842036 * pow(beta, 4) + 108.53608309082300208 * pow(beta, 5)
    + 299.26769419241406575 * pow(beta, 6) + 824.59140866646260666 * pow(beta, 7) + 2256.9464691160956608 * pow(beta, 8)
    + 6174.4168460205032479 * pow(beta, 9) + 16819.879385593690953 * pow(beta, 10) + 45803.222727034040360 * pow(beta, 11)
    + 124363.26835432977776 * pow(beta, 12) + 337573.74963124787949 * pow(beta, 13) + 914347.51881247417342 * pow(beta, 14)
    + 2476042.0363235919004 * pow(beta, 15) + 6694111.3933182426254 * pow(beta, 16) + 18094604.163418613844 * pow(beta, 17)
    + 48847832.893538572297 * pow(beta, 18) + 131848611.02423050678 * pow(beta, 19) + 355511932.47075480765 * pow(beta, 20);

    return x;
}



int main() {
    // HERE IS ONLY USER INPUT //

    //CHOOSE SIZE OF CUBE, AMOUNT OF MONTECARLO STEPS & BIN SIZE FOR JACKKNIFE RESAMPLING
    int size;
    cout << "Enter Lattice Size (int): ";
    cin >> size;

    if (checkCinFail() || size < 0) {
        // Handle the error or exit
        return 1; // Exit the program indicating failure
    }
    // If the input was successful, proceed with the rest of the program
    cout << "Entered Lattice Size: " << size << endl;


    int n_steps;
    cout << "Enter Amount of MC steps(int): ";
    cin >> n_steps;

    if (checkCinFail() || n_steps < 0 || n_steps % 2 != 0) {
        // Handle the error or exit
        return 1; // Exit the program indicating failure
    }

    cout << "Chosen Amount of MC Steps: " << n_steps << endl;

    int N_bin;
    cout << "Enter Amount of Jackknifef Bins (int): ";
    cin >> N_bin;

    if (checkCinFail() || N_bin < 0 || N_bin % 2 != 0) {
        // Handle the error or exit
        return 1; // Exit the program indicating failure
    }

    cout << "Chosen Jackknife Bin Size: " << N_bin << endl;

    int bin_size = n_steps*0.8 / N_bin;

    //TEMPERATURE
    string a;
    cout << "Do you want to study the (3D) Ising Model or Blume-Capel Model at it's T_C (I/BC): ";
    cin >> a;
    double temp;
    double D;

    if (checkCinFail()) {
        // Handle the error or exit
        return 1; // Exit the program indicating failure
    }
    char hightemp;

    if(a == "I" || a == "i") {
        temp = 1/0.22165463;
        D = -1000;
        cout << "Chosen Model: Ising Model (Limit D -> -infinity)" << endl;
    }
    else if(a == "BC" || a == "bc") {
        temp = 1/0.387721735;
        D = 0.655;
        cout << "Chosen Model: Blume-Capel Model with D = 0.655" << endl;
        cout << "Do you want to change the value of D? (Y/N)" << endl;
        char check;
        cin >> check;
        if (check == 'y' || check == 'Y') {
            cout << "Enter the new value for D: " << endl;
            cin >> D;
            cout << "New value of D: " << D << endl;
        }
        cout << "Do you want to compare the value of the susceptibility with the high temperature expansion? (Y/N)" << endl;
        cout << "Notice: This only makes sense if you choose a (inverse) Temperature ß < 0.15 in the next step and D will be set to D=0.641 (Y/N): " << endl;
        cin >> hightemp; if(checkCinFail()) {return 1;}
        if (hightemp == 'Y' || hightemp == 'y') {
            cout << "Update the inverse temperature (beta) to a value < 0.15 (double):  " << endl;
            cin >> temp; if(checkCinFail()) {return 1;}
            temp = 1/temp;
            cout << "The new inverse Temperature is: " << 1/temp << endl;
            D = 0.641;
            cout << "The new value for D is: " << D << endl; 
        }
    }
    
    else{
        cerr << "Invalid input. Program will terminate." << endl;
        return 1; // Indicates an error occurred
    }



    //Periodic Boundaries in xyz-Direction (TRUE/FALSE)
    //Left True statement is for indices at i/j/k = 0, Right True Statement at indices i/j/k = (L-1)

    bool bc_x1, bc_x2, bc_y1, bc_y2, bc_z1, bc_z2;

    char bc;
    cout << "Do you want periodic or open boundary conditions? (P/O): " << endl;
    cin >> bc;
    if (checkCinFail()) {
        // Handle the error or exit
        return 1; // Exit the program indicating failure
    }
    if (bc == 'p' || bc == 'P') {
        bc_x1 = true; bc_x2 = true; bc_y1 = true; bc_y2 = true; bc_z1 = true; bc_z2 = true;
    }
    else if(bc == 'O' || bc == 'o') {
        cout << "I recommend choosing only +z and -z to be open, as this code in it's current form can only calculate the Correlations Surface-Bulk and Surface-Surface if both are turned to Open. Also the average is only taken across these two surfaces." << endl;
        cout << "The Correlations will only be calculated and stored if (+/-x , +/- y, +/- z) = 1 or (+/- x, +/-y) = 1 and (+/- z) = 0)" << endl;
        cout << "You can now turn on/off periodic boundary conditions for any direction you like." << endl;
        cout << "Turn ON = 1; Turn OFF = 0" << endl;
        cout << "+x Direction: " << endl;
        cin >> bc_x1; if (checkCinFail()) {return 1;}
        cout << "-x Direction: " << endl;
        cin >> bc_x2; if (checkCinFail()) {return 1;}
        cout << "+y Direction: " << endl;
        cin >> bc_y1; if (checkCinFail()) {return 1;}
        cout << "-y Direction: " << endl;
        cin >> bc_y2; if (checkCinFail()) {return 1;}
        cout << "+z Direction: " << endl;
        cin >> bc_z1; if (checkCinFail()) {return 1;}
        cout << "-z Direction: " << endl;
        cin >> bc_z2; if (checkCinFail()) {return 1;}
    }
    else {return 1;}
    cout << "The following Directions have periodic boundary conditions (0 = OPEN, 1 = PERIODIC): " << endl;
    cout << "+x: " << bc_x1 << endl;
    cout << "-x: " << bc_x2 << endl;
    cout << "+y: " << bc_y1 << endl;
    cout << "-y: " << bc_y2 << endl;
    cout << "+z: " << bc_z1 << endl;
    cout << "-z: " << bc_z2 << endl;

    pair<const bool,const bool> periodic_boundaries_x (bc_x1,bc_x2);
    pair<const bool,const bool> periodic_boundaries_y (bc_y1,bc_y2);
    pair<const bool,const bool> periodic_boundaries_z (bc_z1,bc_z2);

    //Create Triple of 3 bool pairs, for Object Creation
    triple<pair<const bool,const bool>,pair<const bool,const bool>,pair<const bool,const bool>> Boundaries(periodic_boundaries_x, periodic_boundaries_y, periodic_boundaries_z);

    //Blume-Capel-Parameters
    char extfield;
    double h;
    cout << "Do you want an external magnetic field (Y/N): " << endl;
    cin >> extfield; if (checkCinFail()) {return 1;}
    if (extfield == 'Y' || extfield =='y') {
        cout << "Enter the value for h (double): " << endl;
        cin >> h; if (checkCinFail()) {return 1;}
    }
    else if(extfield == 'N' || extfield =='n') {
        cout << "The external magnetic field is turned off, h = 0; " << endl;
        h = 0.0;
    }
    else{
        return 1;
    }

    char other_observables;
    bool other;
    if (hightemp != 'Y' || hightemp != 'y') {
        cout << "Do you want to calculate other Observables (Magnetization, Susceptibility, Binder Cumulant) (Y/N): " << endl;
        cout << "This is on by default, if you calculate the high-temp expansion of the Susceptibility for the Blume-Capel Model.(Y/N)" << endl;
        cin >> other_observables; if (checkCinFail()) {return 1;}
        if (other_observables == 'y' || other_observables == 'Y') {
            other = true;
        }
        else if (other_observables == 'n' || other_observables == 'N') {
            other = false;
        }
        else{return 1;}
    }
    else {
        other = true;
    }



    //MONTECARLO LOOP, X is used if one would like to calculate the observables for a variety of temperatures,
    //Updating the temperature in each iteration, for example to find the phase transition.
    cout << "Monte Carlo Simulation starting..." << endl;
    //Start of the MonteCarlo Simulation
    for(int x=0; x<1; x++) {
        BlumeCapel q(size, temp, D, h, Boundaries);
        
        //Initial Arrays for Magnetization/Susceptibility Mean
        vector<vector<double>> corr_z(size);
        vector<vector<double>> corr_xy(size);
        vector<vector<double>> corr_bulk(size);

        vector<double> M_vec;
        vector<double> X_vec;
        vector<double> Energy_vec;

        //TRACK TIME
        auto start = chrono::high_resolution_clock::now();

        for (int step=0; step < n_steps; step++) {
            for (int p=0; p<pow(size,3);p++) {
              q.spinFlip();
            }
            for (int l=0; l<size; l++) {
               q.clusterFlip();
            }
            if(step >= (n_steps - n_steps*0.8)) { //THROWS AWAY FIRST 20% OF VALUES
                if(hightemp != 'Y' || hightemp != 'y'){
                    if(periodic_boundaries_z.first == false && periodic_boundaries_z.second == false) {
                        for(int z = 0; z < size - 1; z++) {
                            double sum = 0.0;
                            for(int kx = 0; kx < size; kx++) {
                                for (int ky = 0; ky < size; ky++) {
                                    double v1, v2;
                                    triple<int,int,int> spin_1 (kx,ky,0);
                                    triple<int,int,int> spin_2 (kx,ky,z+1);
                                    triple<int,int,int> spin_3 (kx,ky,size-1);
                                    triple<int,int,int> spin_4 (kx,ky,size-z-2);
                                    v1 = q.tp_correlation(spin_1, spin_2);
                                    v2 = q.tp_correlation(spin_3, spin_4);
                                    sum += (v1+v2)/2;
                                }
                            }
                            sum *= pow(size,-2);
                            corr_z[z].push_back(sum);
                        }
                        for(int y = 0; y < size - 1; y++) {
                            double sum = 0.0;
                            for(int x = 0; x < size; x++) {
                                double v1, v2;
                                triple<int,int,int> spin_1 (x,0,0);
                                triple<int,int,int> spin_2 (x,y + 1,0);
                                triple<int,int,int> spin_3 (x,0,size-1);
                                triple<int,int,int> spin_4 (x,y + 1 ,size-1);
                                v1 = q.tp_correlation(spin_1, spin_2);
                                v2 = q.tp_correlation(spin_3, spin_4);
                                sum += (v1+v2)/2;
                            }
                            sum *= pow(size,-1);
                            corr_xy[y].push_back(sum);
                        }
                    }
                    if(periodic_boundaries_z.first == true && periodic_boundaries_z.second == true) {
                        for(int d = 0; d < size-1; d++) {
                            double sum = 0.0;
                            for(int z = 0; z < size; z++) {
                                for (int y = 0; y < size; y++) {
                                    for (int x = 0; x < size; x++) {
                                        double v1, v2, v3, v4, v5, v6;
                                        triple<int,int,int> spin_1 (x,y,z);
                                        triple<int,int,int> spin_2 (x,y,z+d+1);
                                        v1 = q.tp_correlation(spin_1,spin_2);

                                        triple<int,int,int> spin_3 (x,y,z);
                                        triple<int,int,int> spin_4 (x,y+d+1,z);
                                        v2 = q.tp_correlation(spin_3,spin_4);

                                        triple<int,int,int> spin_5 (x,y,z);
                                        triple<int,int,int> spin_6 (x+d+1,y,z);
                                        v3 = q.tp_correlation(spin_5,spin_6);

                                        triple<int,int,int> spin_7 (x,y,z);
                                        triple<int,int,int> spin_8 (x-d-1,y,z);
                                        v4 = q.tp_correlation(spin_7,spin_8);

                                        triple<int,int,int> spin_9 (x,y,z);
                                        triple<int,int,int> spin_10 (x,y-d-1,z);
                                        v5 = q.tp_correlation(spin_9,spin_10);

                                        triple<int,int,int> spin_11 (x,y,z);
                                        triple<int,int,int> spin_12 (x,y,z-d-1);
                                        v6 = q.tp_correlation(spin_11,spin_12);

                                        sum += (v1+v2+v3+v4+v5+v6)/6;
                                    }
                                }
                            }
                            sum *= pow(size,-3);
                            corr_bulk[d].push_back(sum);
                        }
                    }
                }
                if(other == true) {
                        pair<double,double> c = q.calcMagChi();
                        M_vec.push_back(c.first);
                        X_vec.push_back(c.second);
                        double Energy_val = q.getE()/(pow(size,3));
                        Energy_vec.push_back(Energy_val);
                    }
            }
            if (step % 100 == 0) {
                cout << step << "th step done" << endl;
            }
        }

        if(periodic_boundaries_z.first == false && periodic_boundaries_z.second == false) {

            //PRINT RESULTS IN .TXT FILE
            ofstream c("correlation_surface-bulk_" + to_string(size) + ".txt");
            ofstream c_s("correlation_surface-surface_" + to_string(size) + ".txt");
            c << "z" <<";" << "corr" << ";" << "sig_corr";
            c_s << "y" << ";" << "corr" << ";" << "sig_corr";

            vector<pair<double,double>> corr;
            vector<pair<double,double>> corr_surface;
            for(int z = 0; z < size - 1; z++) {
                pair<double, double> correlation = jack_mean_std(jack_avg(CalcBins(create_bins(corr_z[z], N_bin, bin_size)), 1));
                corr.push_back(correlation);
                pair<double, double> correlation_surface = jack_mean_std(jack_avg(CalcBins(create_bins(corr_xy[z], N_bin, bin_size)), 1));
                corr_surface.push_back(correlation_surface);
            }

            int z=1;
            for (auto i = corr.begin(); i != corr.end(); i++) {
                c << "\n" << z << ";" << (*i).first << ";" << (*i).second;
                z +=1;
            }

            int y=1;
            for(auto j = corr_surface.begin(); j != corr_surface.end(); j++) {
                c_s << "\n" << y << ";" << (*j).first << ";" << (*j).second;
                y += 1;
            }

            c.close();
            c_s.close();
        }

        if(periodic_boundaries_z.first == true && periodic_boundaries_z.second == true) {
            ofstream c_b("correlation_bulk-bulk_2_" + to_string(size) + ".txt");
            c_b << "z';corr;sig_corr";

            vector<pair<double,double>> c_bulk;
            for (int d = 0; d < size - 1; d++) {
                pair<double,double> correlation = jack_mean_std(jack_avg(CalcBins(create_bins(corr_bulk[d], N_bin, bin_size)), 1));
                c_bulk.push_back(correlation);
            }

            int d=1;
            for (auto k = c_bulk.begin(); k != c_bulk.end(); k++) {
                c_b << "\n" << d << ";" << (*k).first << ";" << (*k).second;
                d += 1;
            }
            c_b.close();
        }

        if (other == true) {
            vector<vector<double>> bin_matrix = create_bins(M_vec, N_bin, bin_size);

            
            triple<vector<double>, vector<double>, vector<double>> M = CalcBins_Mag(bin_matrix);
            vector<double> M1 = jack_avg(M.first,1);
            vector<double> M2 = jack_avg(M.second,1);
            vector<double> M4 = jack_avg(M.third,1);

            pair<double,double> U4 = U4_calc(M2, M4);

            // Mean and STD of M and X:
            pair<double,double> M1jack = jack_mean_std(M1);
            pair<double, double> Xjack = jack_mean_std(jack_avg(CalcBins(create_bins(X_vec, N_bin, bin_size)),1));
            pair<double, double> E_f = jack_mean_std(jack_avg(CalcBins(create_bins(Energy_vec, N_bin, bin_size)), 1));

            ofstream f("other_observables.txt");
            f << "U4;U4_err;M;M_err;X;X_err;E;E_err;" << endl;
            f <<  U4.first << ";" << U4.second << ";" << M1jack.first << ";" << M1jack.second << ";" << Xjack.first
            << ";" << Xjack.second << ";" << E_f.first << ";" << E_f.second << endl;
            f.close();

            if(hightemp == 'y' || hightemp == 'Y') {
                double corr_x = expansion(temp);
                cout << "Comparison of the Susceptibility to High-Temp Expansion: " << endl << endl;
                cout << "Calculated Value: " << Xjack.first << "+/-" << Xjack.second << endl;
                cout << "Expansion Value: " << corr_x << endl;
                cout << "Deviation: " << abs(100 - Xjack.first / corr_x * 100) << "%" << endl << "Deviation in Sigma:" << abs(Xjack.first - corr_x)/Xjack.second << endl << endl;
            }
        }



        auto end = chrono::high_resolution_clock::now();
        //CALCULATE SIMULATION TIME:
        auto duration = chrono::duration_cast<chrono::milliseconds>(end-start).count();
        std::cout << "The " << x+1 << "th run took: " << duration/1000.0 << " s to run." << endl << endl;


    }


    return 0;
}