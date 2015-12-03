//
//  main.cpp
//  Grid_diffuse_reflection
//
//  Created by Виктор Ливинюк on 26.11.15.
//  Copyright © 2015 Viktor Liviniuk. All rights reserved.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;
// 2dim Grid (1 space dimention, 1 velocity dimension)
class GridDR {
    int N;          // space grid size
    int N_xi;       // velocity grid size
    double a, b;    // space grid limits
    double xi_cut;  // velocity grid limit (to both sides of 0)
    double xi1;     // sqrt(k * T1 / m)
    double h;       // space grid width of cell
    double dxi;     // velocity grid width of cell
    double dt;      // time step
    
    double C;       // constant in Maxwell distribution
    
    const double k = 1.380648528e-28;       // Boltzmann constant
    const double m = 4.002602 * 1.6605402e-27;
    
    // thermos parameters
    double T1;
    double T2;
    
    // optimization
    double sum2, sum2_;
    
    vector<vector<double> > f;  // function
public:
    GridDR(int const N_, int const N_xi_, double const a_, double const b_, double const T1_, double const T2_) {
        T1 = T1_;
        T2 = T2_;
        N = N_;
        N_xi = N_xi_;
        f = vector<vector<double> >(N + 2, vector<double>(N_xi, 0));
        a = a_;
        b = b_;
        xi1 = sqrt(k * T1 / m);
        xi_cut = 4 * xi1;
        h = (b - a) / N;
        dxi = 2 * xi_cut / N_xi;
        dt = h / xi_cut;            // in this case the Currant condition is satisfied for all velocities
        C = 1;
        
        sum2 = 0;
        for (b = N_xi / 2; b < N_xi; b++)
            sum2 += xi(b) * exp(-m * pow (xi(b), 2) / (2 * k * T1));    // here always xi(b) > 0
        sum2_ = 0;
        for (b = 0; b < N_xi / 2; b++)
            sum2_ += - xi(b) * exp(-m * pow (xi(b), 2) / (2 * k * T2));  // here always xi(b) < 0
        
    }
    
    // get a cell coordinate by index
    // ______________________
    // i actually means (i + 1/2)
    // x(1) == a + h / 2
    // x(N) == b - h / 2
    double x(int const i) { return a + h * (i - 0.5); }
    
    // xi(0)            == -xi_cut + 1/2 * dxi
    // xi(N_xi / 2 - 1) < 0
    // xi(N_xi / 2)     > 0
    // xi(N_xi - 1)     == xi_cut - 1/2 * dxi
    double xi(int const a) { return -xi_cut + dxi * (a + 0.5); }
    
    vector<double> & operator [] (int const i) { return f[i]; }
    
    // n0 - initial concentration
    void set_initial_condition_maxwell(double n0) {
        double c = 0;
        for (int a = 0; a < N_xi; a++)
            c += exp(-m * pow(xi(a), 2) / (2 * k * T2));
        C = n0 / dxi / c;
        for (int i = 1; i <= N; i++)
            for (int a = 0; a < N_xi; a++)
                f[i][a] = exp(-m * pow(xi(a), 2) / (2 * k * T2));
    }
    
    void output(const char* file_name) {
        std::ofstream fout;
        fout.open(file_name);
        for (int i = 1; i <= N; i++) {
            for (int a = 0; a < N_xi; a++) {
                fout << x(i) << ' ' << xi(a) << ' ' << C * f[i][a] << endl;
            }
            fout << endl;
        }
        fout.close();
    }
    
    double approx(double doubled, double negative) {
        double ret = 2 * doubled - negative;
        if (ret < 0)
            return 0;
        return ret;
    }
    
    double df_i(double f_minus, double f, double f_plus) {
        double t1 = f_plus - f;
        double t2 = f - f_minus;
        double a1 = abs(t1);
        double a2 = abs(t2);
        double df = t1 * t2 <= 0 ? 0 : ((a1 < a2 ? a1 : a2) * (t1 > 0 ? 1 : -1));
        return df;
    }
    
    double i_plus_half_minmod_positive(double f_minus, double f, double f_plus, double y) {
        return f + (1 - y) / 2 * df_i(f_minus, f, f_plus);
    };
    
    double i_plus_half_minmod_negative(double f, double f_plus, double f_plus_two, double y) {
        return f_plus - (1 - y) / 2 * df_i(f, f_plus, f_plus_two);
    };
    
    // time++
    void next_step() {
        //virtual nodes
        // starting with v[0] == f[1/2] to v[N] == f[N + 1/2]
        // v[i] == f[i + 1/2]
        vector<vector<double> > v(N + 1, vector<double>(N_xi, 0));
        
        // D
        for (int a = N_xi / 2; a < N_xi; a++) {     // for positive xi(a)
            f[N + 1][a] = approx(f[N][a], f[N - 1][a]);
            v[N][a] = i_plus_half_minmod_positive(f[N - 1][a], f[N][a], f[N + 1][a], xi(a) * dt / h);
        }
        for (int a = 0; a < N_xi / 2; a++) {     // for negative xi(a)
            f[0][a] = approx(f[1][a], f[2][a]);
            v[0][a] = i_plus_half_minmod_negative(f[0][a], f[1][a], f[2][a], xi(a) * dt / h);
        }
        // E
        double sum1 = 0;
        for (b = 0; b < N_xi / 2; b++)
            sum1 += - xi(b) * v[0][b];                                  // here always xi(b) < 0
        for (int a = N_xi / 2; a < N_xi; a++) {
            double e = exp(-m * pow (xi(a), 2) / (2 * k * T1));
            v[0][a] = sum1 * e / sum2;
        }
        
        sum1 = 0;
        for (b = N_xi / 2; b < N_xi; b++)
            sum1 += xi(b) * v[N][b];                                    // here always xi(b) > 0
        for (int a = 0; a < N_xi / 2; a++) {
            double e = exp(-m * pow (xi(a), 2) / (2 * k * T2));
            v[N][a] = sum1 * e / sum2_;
        }
        
        // J
        vector<double> g(N_xi, 0);
        
        sum1 = 0;
        for (b = 0; b < N_xi / 2; b++)
            sum1 += - xi(b) * (f[1][b] - 0.5 * df_i(f[0][b], f[1][b], f[2][b]));
        for (int a = N_xi / 2; a < N_xi; a++) {
            // find g[a]
            double e = exp(-m * pow (xi(a), 2) / (2 * k * T1));
            g[a] = sum1 * e / sum2;
            f[0][a] = approx(g[a], f[1][a]);
            v[1][a] = i_plus_half_minmod_positive(f[0][a], f[1][a], f[2][a], xi(a) * dt / h);
        }
        
        sum1 = 0;
        for (b = N_xi / 2; b < N_xi; b++)
            sum1 += xi(b) * (f[N][b] + 0.5 * df_i(f[N - 1][b], f[N][b], f[N + 1][b]));
        for (int a = 0; a < N_xi / 2; a++) {
            // find g[a]
            double e = exp(-m * pow (xi(a), 2) / (2 * k * T2));
            g[a] = sum1 * e / sum2_;
            f[N + 1][a] = approx(g[a], f[N][a]);
            v[N - 1][a] = i_plus_half_minmod_negative(f[N - 1][a], f[N][a], f[N + 1][a], xi(a) * dt / h);
        }
        
        // we already found elements of v:
        // for xi > 0: v[0][a], v[1][a]     and v[N][a]
        // for xi < 0: v[0][a], v[N - 1][a] and v[N][a]
        // now we seek for all other elements of v
        for (int a = N_xi / 2; a < N_xi; a++) {
            for (int i = 2; i < N; i++) {
                v[i][a] = i_plus_half_minmod_positive(f[i - 1][a], f[i][a], f[i + 1][a], xi(a) * dt / h);
            }
        }
        for (int a = 0; a < N_xi / 2; a++) {
            for (int i = 1; i < N - 1; i++) {
                v[i][a] = i_plus_half_minmod_negative(f[i][a], f[i + 1][a], f[i + 2][a], xi(a) * dt / h);
            }
        }
        
        // now we found all v[][] elements
        // let's find next time-layer
        // auxiliary 2-dim vector
        vector<vector<double> > F(N + 2, vector<double>(N_xi, 0));
        for (int a = 0; a < N_xi; a++) {
            for (int i = 1; i <= N; i++) {
                F[i][a] = f[i][a] - (xi(a) * dt / h) * (v[i][a] - v[i-1][a]);
            }
        }
        f.clear();
        f = F;
        
    }
    
    void play_for(int n) {
        for (int i = 0; i < n; i++)
            next_step();
    }
    
    double temperature(int i) {
        if (i < 0 || i > N)
            return -1;
        double sum1 = 0, sum2 = 0;
        for (int a = 0; a < N_xi; a++) {
            sum1 += pow(xi(a), 2) * f[i][a];
            sum2 += f[i][a];
        }
        if (sum2 < 1e-15)
            return 0;
        return m / k * sum1 / sum2;
    }
};

int main(int argc, const char * argv[]) {
    // thermos parameters
    double T1 = 350;
    double T2 = 200;
    double n0 = 300;
    
    // grid parameters
    int N_ = 20;
    int N_xi_ = 30;
    double a_ = 0;
    double b_ = 10;
    
    // create a thermos
    GridDR grid(N_, N_xi_, a_, b_, T1, T2);
    
    // some actions
    grid.set_initial_condition_maxwell(n0);

//    grid.play_for(1250);
    // output for testing maxwell
    grid.output("/Users/Viktor/Documents/MIPT_base/Grid_diffuse_reflection/Grid_diffuse_reflection/start.txt");
    
    // check correctness
    double T_theor = sqrt(T1 * T2);
    cout << "Theoretical temperature:\t" << T_theor << endl;
    cout << "Grid temperature in   1:\t" << grid.temperature(1) << endl;
    cout << "Grid temperature in N/2:\t" << grid.temperature(N_ / 2) << endl;
    cout << "Grid temperature in   N:\t" << grid.temperature(N_) << endl;

    return 0;
}
