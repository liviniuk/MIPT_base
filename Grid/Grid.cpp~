#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <cmath>

class Grid {
private:
	std::vector< double > data;
	int size;
	double dt;
	double dh;
	double xi;
	void vector_copy(std::vector< double >&, std::vector< double >);
	int steps;
	double i_plus_half_minmod(double, double, double, double);
public:
	Grid();
	Grid(int);
	Grid(int, double, double, double);
	Grid(std::vector< double >, double, double, double);
	void set_size(int);
	void make_pillar(int, int);
	void make_bell(int, int);
	void output(const char*);
	Grid solve(int);
	Grid solve_clear_4point(int);
	Grid solve_unclear_6point(int);
	Grid solve_difference_3layer(int);
	Grid solve_with_minmod_limiter(int);
	Grid exact_solution(int);
	double& operator [] (size_t i) { return data[i]; };
};


Grid::Grid() : dt(1.0), dh(1.0), xi(0.5), size((int) (100)) { data = std::vector< double >(size, 0); }

Grid::Grid(int n) : dt(1.0), dh(1.0), xi(0.5), size(n) { data = std::vector< double >(size, 0); }

//Grid::Grid(int n, double dh_) : dt(1.0), dh(dh_), xi(dh_ / 2.0), size((int) (n / dh_)) { data = std::vector< double >(size, 0); }
Grid::Grid(int n, double dh_, double dt_, double xi_) : dt(dt_), dh(dh_), xi(xi_), size((int) (n / dh_)) { data = std::vector< double >(size, 0); }

Grid::Grid(std::vector< double > vec, double dh_, double dt_, double xi_) : data(vec), dt(dt_), dh(dh_), xi(xi_), size(vec.size()) { 
	data = std::vector< double >(size, 0);
	for (int i = 0; i < size; i++)
		data[i] = vec[i];
}

void Grid::set_size(int num) { 
	size = (int) (num / dh); 
	data.assign(size, 0);
}

void Grid::make_pillar(int a, int b) {
	a = (int) (a / dh);
	b = (int) (b / dh);
	data.assign(size, 0);
	for (int i = a; i < b; i++)
		 data[i] = 1;
}

void Grid::make_bell(int a, int b) {
	a = (int) (a / dh);
	b = (int) (b / dh);
	data.assign(size, 0);
	if ((b - a) % 2 == 0)
		b++;
	int center = (b + a) / 2;
	double c = exp(- (a - center) * (a - center) * dh * dh / 2.0);
	for (int i = a; i < b; i++)
		 data[i] = exp(- (i - center) * (i - center) * dh * dh / 2.0) - c;
}

void Grid::output(const char* file_name) {
	std::ofstream fout;
	fout.open(file_name);
	for (int i = 0; i < size; i++)
		fout << i * dh << ' ' << data[i] << std::endl;
	fout.close();
}

void Grid::vector_copy(std::vector< double >& vec, std::vector< double > source) {
	vec.clear();
	for(int i = 0; i < size; i++)
		vec.push_back(source[i]);
}

Grid Grid::solve(int time_steps) {
	std::vector< double > data(this->data);	
	std::vector< double > auxiliary;
	int steps = time_steps / dh;
	for (int i = 0; i < steps; i++) {
		vector_copy(auxiliary, data);
		for (int j = 1; j < size; j++) {
			data[j] = auxiliary[j] - xi * dt * (auxiliary[j] - auxiliary[j - 1]) / dh;
		}
	}
	return Grid(data, dh, dt, xi);
}

Grid Grid::solve_clear_4point(int time_steps) {
	std::vector< double > data(this->data);	
	std::vector< double > auxiliary;
	double y = xi * dt / dh;
	double t1 = y / 2 * (1 + y);
	double t2 = (1 - y * y);
	double t3 = y / 2 * (1 - y);
	int steps = time_steps / dh;
	for (int i = 0; i < steps; i++) {
		vector_copy(auxiliary, data);
		for (int j = 1; j < size; j++) {
			data[j] = t1 * auxiliary[j - 1] + t2 * auxiliary[j] - t3 * auxiliary[j + 1];
		}
	}
	return Grid(data, dh, dt, xi);
}

Grid Grid::solve_unclear_6point(int time_steps) {
	std::vector< double > data(this->data);	
	std::vector< double > auxiliary;
	std::vector< double > alpha;
	std::vector< double > beta;
	
	double y = dt * xi / dh;		
	int steps = time_steps / dh;
	for (int j = 0; j < steps; j++) {
		vector_copy(auxiliary, data);
		alpha.assign(size, 0);
		beta.assign(size, 0);
		data[size - 1] = data[0] = 0;
		for (int i = 1; i < size - 1; i++) {
			double D_i = auxiliary[i] - y / 4 * (auxiliary[i + 1] - auxiliary[i - 1]);
			alpha[i + 1] = 1 / (alpha[i] - 4 / y);
			beta[i + 1] = - alpha[i + 1] * (beta[i] + 4 * D_i / y);
		}
		for (int i = size - 2; i >= 1 ; i--) 
			data[i] = alpha[i + 1] * data[i + 1] + beta[i + 1];
	}
	return Grid(data, dh, dt, xi);
}

Grid Grid::solve_difference_3layer(int time_steps) {
	std::vector< double > data(this->data);	
	std::vector< double > auxiliary1;
	std::vector< double > auxiliary2(size, 0);
	for (int i = 0; i < size - 1; i++)
		auxiliary2[i + 1] = data[i];

	double y = dt * xi / dh;		
	int steps = time_steps / dh;
	for (int j = 1; j < steps; j++) { // from 1 becouse we start with the 3rd layer, not with the 2nd one
		vector_copy(auxiliary1, auxiliary2);
		vector_copy(auxiliary2, data);
		data[0] = data[size - 2] = 0;
		for (int i = 1; i < size - 1; i++) 
			data[i] = auxiliary1[i] + y * (auxiliary2[i - 1] - auxiliary2[i + 1]);
	}
	return Grid(data, dh, dt, xi);
}

double Grid::i_plus_half_minmod(double f_minus, double f, double f_plus, double y) {
	double t1 = f_plus - f;
	double t2 = f - f_minus;
	double df = t1 * t2 <= 0 ? 0 : ((t1 < t2 ? t1 : t2) * (t1 > 0 ? 1 : -1));
	return f + (1 - y) / 2 * df;
};

Grid Grid::solve_with_minmod_limiter(int time_steps) {
	std::vector< double > data(this->data);	
	std::vector< double > auxiliary;
	std::vector< double > halfs; // halfs[i] == f[i + 1/2]

	double y = dt * xi / dh;		
	int steps = time_steps / dh;
	for (int j = 0; j < steps; j++) {
		vector_copy(auxiliary, data);
		//find halfs
		halfs.clear();
		halfs.push_back(0); //halfs[0] = f[1/2]
		for (int i = 1; i < size - 2; i++) {
			double f_i_plus_half = i_plus_half_minmod(auxiliary[i - 1], auxiliary[i], auxiliary[i + 1], y);
			halfs.push_back(f_i_plus_half);
		}
		halfs.push_back(0); // halfs[size - 2] == f[size - 3/2] = 0
		data[0] = data[size - 1] = 0;
		for (int i = 1; i < size - 1; i++)
			data[i] = auxiliary[i] + y * (halfs[i - 1] - halfs[i]);
	}
	return Grid(data, dh, dt, xi);
}

Grid Grid::solve_diffusion(int time_steps) {
	int N = 20;
	double a = 0, b = 5;
	double h = (b - a) / N;
	int N_xi = 20;
	double k = 1.3806485e-23, T = 300, m = 4.002602 * 1.6605402e-27;
	double xi_kv = sqrt(k * T / m);
	double xi_cut = 5 * xi_kv;
	double dxi = 2 * xi_cut / N_xi;
	//??double dt = h / xi_cut;
	std::vector< double > w_h, w_xi;
	double temp = a;
	for (int i = 0; i < N_xi; i++)
		w_h.push_back(h * (i - 0.5));
	for (int i = 0; i < N_xi; i++)
		w_xi.push_back(-xi_cut + dxi * (i - 0.5));
	std::vector< std::vector< double > > f(N, vector < double >(N_xi, 0));
	// Make a bell
	double p = (int) (5 / dh);
	double q = (int) (8 / dh);
	if ((q - p) % 2 == 0)
		q++;
	int center = (p + q) / 2;
	double c = exp(- (p - center) * (p - center) * dh * dh / 2.0);
	for (int i = p; i < q; i++)
		 f[0][i] = exp(- (i - center) * (i - center) * dh * dh / 2.0) - c;
	
	
	return f;
}

Grid Grid::exact_solution(int time_steps) {
	std::vector<double> sol(size, 0);
	int shift = (int) (dt * xi / dh * time_steps / dh);
	for (int i = 0; i < size - shift; i++)
		sol[i + shift] = data[i];
	return Grid(sol, dh, dt, xi);
};

double error(Grid a, Grid b, int size, double dh) {
	double err = 0;
	for (int i = 0; i < size; i++)
		err += fabs(a[i] - b[i]);
	return err * dh;
}

int main(int argc, char* argv[]) { // [time_steps, grid_size, [a, b]]
	int grid_size = 20;
	int p = 3;
	int q = 8;
	double dh = 0.2;
	double dt = 0.2;
	double xi = dh / 2.0 / dt; // y == xi * dt / dh = 0.5
	int shift = 16;
	int time_steps = shift;
	Grid a(grid_size, dh, dt, xi);
	a.make_bell(p, q);
	a.output("start.txt");
	a.exact_solution(time_steps).output("exact.txt");
	a.solve(time_steps).output("output.txt");
	a.solve_clear_4point(time_steps).output("output_clear_4point.txt");
	a.solve_unclear_6point(time_steps).output("output_unclear_6point.txt");
	a.solve_difference_3layer(time_steps).output("output_difference_3layer.txt");
	a.solve_with_minmod_limiter(time_steps).output("output_minmod_limiter.txt");
	//find error
	std::cout << "error: " << error(a.solve_with_minmod_limiter(time_steps), a.exact_solution(time_steps), grid_size / dh, dh) << std::endl;
	return 0;
}



