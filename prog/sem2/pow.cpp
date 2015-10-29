#include<iostream>
#include<string>

double power(double a, int n) {
	if (a == 0)
		return 0;
	double answer = 1;
	double a_2_k = a;
	for (int k = 0; n != 0; k++) {
		answer *= (n % 2 == 1) ? a_2_k : 1;
		n /= 2;
		a_2_k *= a_2_k;
	}
	return answer;
}

int main(int argc, char* argv[]) {
	std::cout << "Input number and power here, with a space between: ";
	double a;
	int n;
	std::cin >> a >> n;
	std::cout << a << " ^ " << n << " = " << power(a, n) << std::endl; 
	return 0;
}
