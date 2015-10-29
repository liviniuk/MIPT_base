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
	void vector_from_data(std::vector< double >&);
	int steps;
	double i_plus_half_minmod(double, double, double, double);
public:
	Grid();
	Grid(int);
	Grid(int, double);
	Grid(std::vector< double >, double);
	double get(int);
	void set_size(int);
	void make_pillar(int, int);
	void make_bell(int, int);
	void output(const char*);
	void solve(int);
	void solve_clear_4point(int);
	void solve_unclear_6point(int);
	void solve_difference_3layer(int);
	void solve_with_minmod_limiter(int);
	Grid exact_solution(int);
	double error();
};


Grid::Grid() : dt(1.0), dh(1.0), xi(0.5), size((int) (100)) { data = std::vector< double >(size, 0); }

Grid::Grid(int n) : dt(1.0), dh(1.0), xi(0.5), size(n) { data = std::vector< double >(size, 0); }

Grid::Grid(int n, double dh_) : dt(1.0), dh(dh_), xi(dh_ / 2.0), size((int) (n / dh_)) { data = std::vector< double >(size, 0); }

Grid::Grid(std::vector< double > vec, double dh_) : data(vec), dt(1.0), dh(dh_), xi(dh_ / 2.0), size(vec.size()) { 
	data = std::vector< double >(size, 0);
	for (int i = 0; i < size; i++)
		data[i] = vec[i];
}

double Grid::get(int i) { return data[i]; }

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
	double e_a = exp(- (a - center) / 2.0 * (a - center));
	for (int i = a; i < b; i++)
		 data[i] = (exp(- (i - center) * (i - center) / 2.0) - e_a);
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


void Grid::vector_from_data(std::vector< double >& vec) {
	vec.clear();
	for(int i = 0; i < size; i++)
		vec.push_back(data[i]);
}

void Grid::solve(int time_steps) {
	std::vector< double > auxiliary;
	int steps = time_steps / dh;
	for (int i = 0; i < steps; i++) {
		vector_from_data(auxiliary);
		for (int j = 1; j < size; j++) {
			data[j] = auxiliary[j] - xi * dt * (auxiliary[j] - auxiliary[j - 1]) / dh;
		}
	}
	steps = time_steps;
}

void Grid::solve_clear_4point(int time_steps) {
	std::vector< double > auxiliary;
	double y = xi * dt / dh;
	double t1 = y / 2 * (1 + y);
	double t2 = (1 - y * y);
	double t3 = y / 2 * (1 - y);
	int steps = time_steps / dh;
	for (int i = 0; i < steps; i++) {
		vector_from_data(auxiliary);
		for (int j = 1; j < size; j++) {
			data[j] = t1 * auxiliary[j - 1] + t2 * auxiliary[j] - t3 * auxiliary[j + 1];
		}
	}
	steps = time_steps;
}

void Grid::solve_unclear_6point(int time_steps) {
	std::vector< double > auxiliary;
	std::vector< double > alpha;
	std::vector< double > beta;
	
	double y = dt * xi / dh;		
	int steps = time_steps / dh;
	for (int j = 0; j < steps; j++) {
		vector_from_data(auxiliary);
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
	steps = time_steps;
}

void Grid::solve_difference_3layer(int time_steps) {
	std::vector< double > auxiliary1;
	Grid a = Grid(data, dh);
	std::vector< double > auxiliary2 = a.exact_solution(1).data;
	
	double y = dt * xi / dh;		
	int steps = time_steps / dh;
	for (int j = 1; j < steps; j++) { // from 1 becouse we start with the 3rd layer, not with the 2nd one
		vector_copy(auxiliary1, auxiliary2);
		vector_from_data(auxiliary2);
		data[0] = data[size - 2] = 0;
		for (int i = 1; i < size - 1; i++) 
			data[i] = auxiliary1[i] + y * (auxiliary2[i - 1] - auxiliary2[i + 1]);
	}
	steps = time_steps;
}

double Grid::i_plus_half_minmod(double f_minus, double f, double f_plus, double y) {
	double t1 = f_plus - f;
	double t2 = f - f_minus;
	double df = t1 * t2 <= 0 ? 0 : ((t1 < t2 ? t1 : t2) * (t1 > 0 ? 1 : -1));
	return f + (1 - y) / 2 * df;
};

void Grid::solve_with_minmod_limiter(int time_steps) {
	std::vector< double > auxiliary;
	std::vector< double > halfs; // halfs[i] == f[i + 1/2]

	double y = dt * xi / dh;		
	int steps = time_steps / dh;
	for (int j = 0; j < steps; j++) {
		vector_from_data(auxiliary);
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
	steps = time_steps;
}

Grid Grid::exact_solution(int time_steps) {
	std::vector<double> sol(size, 0);
	int shift = (int) (dt * xi / dh * time_steps / dh);
	for (int i = 0; i < size - shift; i++)
		sol[i + shift] = data[i];
	return Grid(sol, dh);
};

double Grid::error() {
	Grid exact = exact_solution(steps);
	double err = 0;
	for (int i = 0; i < size; i++)
		err += fabs(data[i] - exact.get(i));
	return err;
}

int main() {
	int time_steps = 16;
	int grid_size = 20;
	Grid a(grid_size, 0.05);
	a.make_pillar(3, 8);
	a.output("start.txt");
	a.exact_solution(time_steps).output("exact.txt");
	a.solve(time_steps);
	a.output("output.txt");
	a.make_pillar(3, 8);
	a.solve_clear_4point(time_steps);
	a.output("output_clear_4point.txt");
	a.make_pillar(3, 8);
	a.solve_unclear_6point(time_steps);
	a.output("output_unclear_6point.txt");
	a.make_pillar(3, 8);
	a.solve_difference_3layer(time_steps);
	a.output("output_difference_3layer.txt");
	a.make_pillar(3, 8);
	a.solve_with_minmod_limiter(time_steps);
	a.output("output_minmod_limiter.txt");
	//std::cout << a.error() << std::endl;
	return 0;
}



