#include <cmath>
#include <fstream>
#include <chrono>
#include "blume.hpp"
#include "statistics.hpp"
using namespace std;


int main() {
    //CHOOSE SIZE OF CUBE, AMOUNT OF MONTECARLO STEPS & BIN SIZE FOR JACKKNIFE RESAMPLING
    int size = 64;
    int n_steps = 10000;
    int N_bin = 100;
    int bin_size = n_steps*0.8 / N_bin;

    //TEMPERATURE
    double temp = 1/0.387721735;

    //Periodic Boundaries in xyz-Direction (TRUE/FALSE)
    //Left True statement is for indices at i/j/k = 0, Right True Statement at indices i/j/k = (L-1)
    pair<const bool,const bool> periodic_boundaries_x (true,true);
    pair<const bool,const bool> periodic_boundaries_y (true,true);
    pair<const bool,const bool> periodic_boundaries_z (true,true);

    //Create Triple of 3 bool pairs, for Object Creation
    triple<pair<const bool,const bool>,pair<const bool,const bool>,pair<const bool,const bool>> Boundaries(periodic_boundaries_x, periodic_boundaries_y, periodic_boundaries_z);

    //Blume-Capel-Parameters
    double D = 0.655;
    double h = 0.0;

    ofstream raw_data_bulk ("NewCalc/raw_data_bulk" + to_string(size) + ".txt");
    raw_data_bulk << "z';corr";
    
    ofstream raw_data_surface_surface("NewCalc/raw_data_surface-surface" + to_string(size) + ".txt");
    raw_data_surface_surface << "y" <<";" << "corr";

    ofstream raw_data_surface_bulk("NewCalc/raw_data_surface-bulk" + to_string(size) + ".txt");
    raw_data_surface_bulk << "z" <<";" << "corr";
    



    //MONTECARLO LOOP
    for(int x=0; x<1; x++) {
        BlumeCapel q(size, temp, D, h, Boundaries);
        
        //Initial Arrays for Magnetization/Susceptibility Mean
        vector<vector<double>> corr_z(size);
        vector<vector<double>> corr_xy(size);

        vector<vector<double>> corr_bulk(size);

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
                if(periodic_boundaries_z.first == false && periodic_boundaries_z.second == false) {
                    for(int z = 0; z < size - 1; z++) {
                        double sum = 0.0;
                        for(int kx = 0; kx < size; kx++) {
                            for (int ky = 0; ky < size; ky++) {
                                triple<int,int,int> spin_1 (kx,ky,0);
                                triple<int,int,int> spin_2 (kx,ky,z+1);
                                triple<int,int,int> spin_3 (kx,ky,size-1);
                                triple<int,int,int> spin_4 (kx,ky,size-z-2);
                                sum += q.tp_correlation(spin_1, spin_2);
                                sum += q.tp_correlation(spin_3, spin_4);
                            }
                        }
                        sum *= 1/2*pow(size,-2);
                        corr_z[z].push_back(sum);
                        raw_data_surface_bulk << "\n" <<  z + 1 << ";" << sum;
                    }
                    for(int y = 0; y < size - 1; y++) {
                        double sum = 0.0;
                        for(int x = 0; x < size; x++) {
                            triple<int,int,int> spin_1 (x,0,0);
                            triple<int,int,int> spin_2 (x,y + 1,0);
                            triple<int,int,int> spin_3 (x,0,size-1);
                            triple<int,int,int> spin_4 (x,y + 1 ,size-1);
                            sum += q.tp_correlation(spin_1,spin_2);
                            sum += q.tp_correlation(spin_3,spin_4);
                        }
                        sum *= 1/2*pow(size,-1);
                        corr_xy[y].push_back(sum);
                        raw_data_surface_surface << "\n" << y + 1 << ";" << sum;
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
                        raw_data_bulk << "\n" << d+1 << ";" << sum;
                    }
                }
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


        auto end = chrono::high_resolution_clock::now();
        //CALCULATE SIMULATION TIME:
        auto duration = chrono::duration_cast<chrono::milliseconds>(end-start).count();
        std::cout << "The " << x+1 << "th run took: " << duration/1000.0 << " s to run." << endl << endl;
    }
    raw_data_bulk.close();
    raw_data_surface_bulk.close();
    raw_data_surface_surface.close();


    return 0;
}