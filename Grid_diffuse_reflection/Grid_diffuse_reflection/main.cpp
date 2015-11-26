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

using namespace std;

// 2dim Grid (1 space dimention, 1 velocity dimension)
class GridDR {
    int N;          // space grid size
    int N_xi;       // velocity grid size
    double a, b;    // space grid limits
    double xi_cut;  // velocity grid limit (to both sides of 0)
    double h;       // space grid width of cell
    double dxi;     // velocity grid width of cell
    double dt;      // time step
    
    const double k = 1.380648528e-28;       // Boltzmann constant
    const double m = 4.002602 * 1.6605402e-27;
    
    // thermos parameters
    double T1;
    double T2;
    
    vector<vector<double> > f;  // function
public:
    GridDR(int const N_, int const N_xi_, double const a_, double const b_, double const xi_cut_) {
        N = N_;
        N_xi = N_xi_;
        f = vector<vector<double> >(N + 2, vector<double>(N_xi, 0));
        a = a_;
        b = b_;
        xi_cut = xi_cut_;
        h = (b - a) / N;
        dxi = 2 * xi_cut / N_xi;
        dt = h / xi_cut;            // in this case the Currant condition is satisfied for all velocities
    }
    
    void set_thermos(double const T1, double const T2) {
        this -> T1 = T1;
        this -> T2 = T2;
    }
    
    // get a cell coordinate by index
    // ______________________
    // i actually means (i + 1/2)
    // x(1) == a + h / 2
    // x(N) == b - h / 2
    int x(int const i) { return a + h * (i - 1/2); }
    // xi(0)        == -xi_cut + 1/2 * dxi
    // xi(N_xi - 1) == xi_cut - 1/2 * dxi
    int xi(int const a) { return -xi_cut + dxi * (a + 1/2); }
    
    // get a virtual cell space coordinate by index
    // y(0) == a
    // y(N) == b
    int y(int const i) { return a + h * i; }
    
    vector<double> & operator [] (int const i) { return f[i]; }
    
    // n0 - initial concentration
    void set_initial_condition_maxwell(double n0) {
        double c = 0;
        for (int a = 0; a < N_xi; a++)
            c += exp(-m * pow(xi(a), 2) / (2 * k * T2));
        c = n0 / dxi / c;
        for (int i = 1; i <= N; i++)
            for (int a = 0; a < N_xi; a++)
                f[i][a] = c * exp(-m * pow(xi(a), 2) / (2 * k * T2));
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
        double df = t1 * t2 <= 0 ? 0 : ((t1 < t2 ? t1 : t2) * (t1 > 0 ? 1 : -1));
        return df;
    }
    
    double i_plus_half_minmod_positive(double f_minus, double f, double f_plus, double y) {
        return f + (1 - y) / 2 * df_i(f_minus, f, f_plus);
    };
    
    double i_plus_half_minmod_negative(double f, double f_plus, double f_plus_two, double y) {
        return f - (1 - y) / 2 * df_i(f, f_plus, f_plus_two);
    };
    
    // time++
    void next_step() {
        //virtual nodes
        // starting with v[0] == f[1/2] to v[N] == f[N + 1/2]
        // v[i] == f[i + 1/2]
        vector<vector<double> > v(N + 1, vector<double>(N_xi, 0));
        
        // D
        for (int a = N_xi / 2; a < N_xi; a++) {
            f[N + 1][a] = approx(f[N][a], f[N - 1][a]);
            v[N][a] = i_plus_half_minmod_positive(f[N - 1][a], f[N][a], f[N + 1][a], xi(a) * dt / h);
        }
        for (int a = 0; a < N_xi / 2; a++) {
            f[0][a] = approx(f[1][a], f[2][a]);
            v[0][a] = i_plus_half_minmod_negative(f[0][a], f[1][a], f[2][a], xi(a) * dt / h);
        }
        // E
        for (int a = N_xi / 2; a < N_xi; a++) {
            double e = exp(-m * pow (xi(a), 2) / (2 * k * T1));
            double sum1 = 0, sum2 = 0;
            for (b = 0; b < N_xi / 2; b++)
                sum1 += - xi(b) * v[0][b];
            for (b = N_xi / 2; b < N_xi; b++)
                sum2 += xi(b) * exp(-m * pow (xi(b), 2) / (2 * k * T1));
            v[0][a] = sum1 * e / sum2;
        }
        for (int a = 0; a < N_xi / 2; a++) {
            double e = exp(-m * pow (xi(a), 2) / (2 * k * T2));
            double sum1 = 0, sum2 = 0;
            for (b = N_xi / 2; b < N_xi; b++)
                sum1 += xi(b) * v[N][b];
            for (b = 0; b < N_xi / 2; b++)
                sum2 += - xi(b) * exp(-m * pow (xi(b), 2) / (2 * k * T2));
            v[N][a] = sum1 * e / sum2;
        }
        
        // J
        vector<double> g(N_xi, 0);
        for (int a = N_xi / 2; a < N_xi; a++) {
            // find g[a]
            double e = exp(-m * pow (xi(a), 2) / (2 * k * T2));
            double sum1 = 0, sum2 = 0;
            for (b = N_xi / 2; b < N_xi; b++)
                sum1 += xi(b) * (f[N][b] + 0.5 * df_i(f[N - 1][b], f[N][b], f[N + 1][b]));
            for (b = 0; b < N_xi / 2; b++)
                sum2 += - xi(b) * exp(-m * pow (xi(b), 2) / (2 * k * T2));
            g[a] = sum1 * e / sum2;
            // find f[0][a]
            f[0][a] = 2 * g[a] - f[1][a];
            if (f[0][a] < 0)
                f[0][a] = 0;
            
            v[1][a] = i_plus_half_minmod_positive(f[0][a], f[1][a], f[2][a], xi(a) * dt / h);
        }
        for (int a = 0; a < N_xi / 2; a++) {
            // find g[a]
            double e = exp(-m * pow (xi(a), 2) / (2 * k * T1));
            double sum1 = 0, sum2 = 0;
            for (b = 0; b < N_xi / 2; b++)
                sum1 += - xi(b) * (f[1][b] - 0.5 * df_i(f[0][b], f[1][b], f[2][b]));
            for (b = N_xi / 2; b < N_xi; b++)
                sum2 += xi(b) * exp(-m * pow (xi(b), 2) / (2 * k * T1));
            g[a] = sum1 * e / sum2;
            // find F[N + 1][a]
            f[N + 1][a] = 2 * g[a] - f[N][a];
            if (f[N + 1][a] < 0)
                f[N + 1][a] = 0;
            
            v[N - 1][a] = i_plus_half_minmod_negative(f[N - 1][a], f[N][a], f[N + 1][a], xi(a) * dt / h);
        }
        
        // we already found elements of v:
        // for xi > 0: v[0][a], v[1][a]     and v[N][a]
        // for xi < 0: v[0][a], v[N - 1][a] and v[N][a]
        // now we seek for all other elements of v
        for (int a = N_xi / 2; a < N_xi; a++) {
            for (int i = 2; i < N; i++) {
                v[i][a] = i_plus_half_minmod_positive(v[i - 1][a], v[i][a], v[i + 1][a], xi(a) * dt / h);
            }
        }
        for (int a = 0; a < N_xi / 2; a++) {
            for (int i = 1; i < N - 1; i++) {
                v[i][a] = i_plus_half_minmod_negative(v[i][a], v[i + 1][a], v[i + 2][a], xi(a) * dt / h);
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
};

int main(int argc, const char * argv[]) {
    // grid parameters
    int N_ = 500;
    int N_xi_ = 500;
    double a_ = 0;
    double b_ = 20;
    double xi_cut_ = 2000;
    
    // thermos parameters
    double T1 = 200;
    double T2 = 300;
    double n0 = 30;
    
    // create a thermos
    GridDR grid(N_, N_xi_, a_, b_, xi_cut_);
    grid.set_thermos(T1, T2);
    
    // some actions
    grid.set_initial_condition_maxwell(n0);

    std::cout << grid[300][300] << endl;
    return 0;
}
