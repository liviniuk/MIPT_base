#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <utility>

class Matrix {
private:
	std::vector< std::vector<double> > data;
	int dim;
public:
	void ran(int);
	int get_dim() const;
	double get(int, int) const;
	void set(int, int, double);
	void operator += (const Matrix &);
	void operator -= (const Matrix &);
	void operator *= (const Matrix &);
	void transpose();
	Matrix();
	Matrix(int);
	Matrix(std::vector<double>);
	Matrix(bool, int);
	void print();
	std::pair<Matrix, Matrix> LU_factorization();
	double determinant();
	Matrix inverse();
};

Matrix::Matrix() : dim(0) { }

Matrix::Matrix(int n) : dim(n), data(n, std::vector<double>(n, 0)) { } 

Matrix::Matrix(std::vector<double> a) : dim(a.size()), data(a.size(), std::vector<double>(a.size(), 0)) {
	for(int i = 0; i < dim; i++)
		data[i][i] = a[i];
}

Matrix::Matrix(bool b, int n) : Matrix(std::vector<double>(n,1)) { }

int Matrix::get_dim() const{ return dim; }

void Matrix::transpose() {
	double buf;
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			buf = data[i][j];
			data[i][j] = data[j][i];
			data[j][i] = buf;
		}
	}
}

void Matrix::operator += (const Matrix &b) {
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			data[i][j] += b.get(i,j);
		}
	}
}

void Matrix::operator -= (const Matrix &b) {
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			data[i][j] -= b.get(i,j);
		}
	}
}
void Matrix::operator *= (const Matrix &b) {
	std::vector<std::vector<double> > buf(dim, std::vector<double>(dim,0));
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			for(int k = 0; k < dim; k++) {
				buf[i][j] += data[i][k] * b.get(k,j);
			}
		}
	}
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			data[i][j] = buf[i][j];
		}
	}
}

double Matrix::get(int i, int j) const { return data[i][j]; }

void Matrix::set(int i, int j, double a) { data[i][j] = a; }

void Matrix::print() {
	for(int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (data[i][j] < 1e-15)
				std::cout << 0 << '\t';			
			else
				std::cout << data[i][j] << '\t';
		}
		std::cout << '\n';
	}
}

void Matrix::ran(int m) {
	srand(time(0));
	for (int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			data[i][j] = (double)(rand() % m); 
		}
	}
}

std::pair<Matrix, Matrix> Matrix::LU_factorization() {
	Matrix L(dim), U(dim);
	for (int j = 0; j < dim; j++) {
		L.set(j, j, 1);
		U.set(0, j, data[0][j]);
	}
	double u_0_0 = U.get(0, 0);
	for (int j = 1; j < dim; j++) {
		L.set(j, 0, data[j][0] / u_0_0);
	}
	for (int i = 1; i < dim; i++) {
		for (int j = i; j < dim; j++) {
			double sum = 0;
			for (int k = 0; k < i;  k++)
				sum += L.get(i, k) * U.get(k, j);
			U.set(i, j, data[i][j] - sum);
		}
		double u_i_i = U.get(i, i);
		for (int j = i + 1; j < dim; j++) {
			double sum = 0;
			for (int k = 0; k < i;  k++)
				sum += L.get(j, k) * U.get(k, i);
			L.set(j, i, (data[j][i] - sum) / u_i_i);
		}
	}
	return std::pair<Matrix, Matrix>(L, U);
}

double Matrix::determinant() {
	double detU = 0;
	std::pair< Matrix, Matrix > p = (*this).LU_factorization();
	for (int i = 0; i < dim; i++) {
		detU += p.second.get(i, i);
	}
	return detU; // detL = 1 (L[i][i] = 1 for all i)
}

Matrix Matrix::inverse() {
	Matrix m(dim);
	std::pair< Matrix, Matrix > p = (*this).LU_factorization();
	for (int i = dim - 1; i >= 0; i--) {
		for (int j = dim - 1; j >= 0; j--) {
			if (i == j == dim - 1)
				m.set(dim - 1, dim - 1, 1 / p.second.get(dim - 1, dim - 1));
			else {
				double sum = 0;
				if (i == j) {
					for (int k = j + 1; k < dim; k++)
						sum += p.second.get(j, k) * m.get(k, j);
					m.set(j, j, (1 - sum) / p.second.get(j, j));
				} else if (i < j) {
					for (int k = i + 1; k < dim; k++)
						sum += p.second.get(i, k) * m.get(k, j);
					m.set(i, j, -sum / p.second.get(i, i));
				} else if (i > j) {
					for (int k = j + 1; k < dim; k++)
						sum += p.first.get(k, j) * m.get(i, k);
					m.set(i, j, -sum);
				}
			}
		}
	}
	return m;
}


int main(int argc, char **argv){
	Matrix a(5);
	std::cout << "\t\tA\n";
	a.print();
	Matrix b(5);
	b.ran(9);
	b += Matrix(std::vector<double>(5, 10));
	std::cout << "\t\tB\n";
	b.print();
	std::pair<Matrix, Matrix> p = b.LU_factorization();
	std::cout << "\t\tL\n";
	p.first.print();
	std::cout << "\t\tU\n";
	p.second.print();
	std::cout << "\t\tLU\n";
	a += p.first;
	a *= p.second;
	a.print();
	std::cout << "\t\tB^-1\n";
	a = b.inverse();
	a.print();
	std::cout << "\t\tB^-1 * B\n";
	a *= b;
	a.print();
	return 0;
}
