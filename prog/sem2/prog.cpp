#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>

bool myCmp(std::string s1, std::string s2) {
	int len1 = s1.length(), len2 = s2.length();
	if (len1 != len2)
		return len1 > len2 ? true : false;
	for (int i = 0; i < len1; i++) {
		if (s1[i] > s2[i])
			return true;
	}
	return false;
}

int main(int argc, char* argv[]) {
	std::vector<std::string> vec;
	std::ifstream file;
	file.open("input.txt");
	if (file != NULL) {
		std::string line;
		while (file >> line) {
			vec.push_back(line);
		}
		std::sort(vec.begin(), vec.end(), myCmp);
		std::ofstream file2;
		file2.open("output.txt");
		if (file2 != NULL) {
			int len = vec.size();
			for (int i = 0; i < len; i++) 	{
				file2 << vec[i] << std::endl;
			}
			file2.close();
		}	
		file.close();
	}
	return 0;
}
