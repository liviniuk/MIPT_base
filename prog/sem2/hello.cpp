#include<iostream>
#include<string>
#include<vector>

int main(int argc, char* argv[]) {
	std::string text("Hello Word!");
	int size = text.size();
	std::vector<char> a(size);
	for (auto i = 0; i < size; i++)
		a[i] = text[i];
	for (auto i = 0; i < size; i++)
		std::cout << a[i];
	std::cout << std::endl;
	return 0;
}
