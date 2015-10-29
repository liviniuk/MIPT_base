#include<iostream>
#include<fstream>
#include<string>
int abs(int a) { return a >= 0 ? a : -a; }

short abs(short a) { return a < 0 ? -a : a; } 

bool polindrom(std::string line) {
	int len = line.length();
	for (int i = 0; i < len / 2; i++) {
		if (line[i] != line [len - i - 1])
			return false;
	}
	return true;
}

int main(int argc, char* argv[]) {
	std::ifstream file;
	file.open("input.txt");
	if (file != NULL) {
		std::string line;
		while (file >> line) {
			if (polindrom(line)) {
				std::cout << "\"" << line << "\"" << " is a polindrom";
			} else {
				std::cout << "\"" << line << "\"" << " is not a polindrom";
			}	
			std::cout << std::endl;
		}
	}
	return 0;
}
