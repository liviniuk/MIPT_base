#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>

template <typename T, typename S>
void copy(T &a, S &b, size_t n) {
	for (size_t i = 0; i < n; i++)
		b[i] = a[i];
}

int  main(int argc, char* argvp[]) {
	const size_t n = 10;
	int a[n];
	for (size_t i = 0; i < n; i++)
		a[i] = i;
	std::vector<int> b(n);
	copy(a, b, n);
	for (size_t i = 0; i < n; i++)
		std::cout << b[i] << std::endl;
}
