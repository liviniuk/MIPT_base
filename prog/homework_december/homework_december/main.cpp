//
//  main.cpp
//  homework_december
//
//  Created by Виктор Ливинюк on 26.11.15.
//  Copyright © 2015 Viktor Liviniuk. All rights reserved.
//

#include <iostream>
#include <map>
#include <fstream>
#include <string>

template<int n>
struct Fibonacci {
    static const int f = Fibonacci<n - 1>::f + Fibonacci<n - 2>::f;
};

template<>
struct Fibonacci<1> {
    static const int f = 1;
};

template<>
struct Fibonacci<2> {
    static const int f = 1;
};

std::map<std::string, int> mapper(std::string filename) {
    std::map<std::string, int> lib;
    std::ifstream file(filename);
    std::string s = "";
    while (std::getline(file, s) ) {
        std::string word;
        for (int i = 0; i < s.length(); i++) {
            if (s[i] == ' ' || s[i] == ',' || s[i] == '.' || s[i] == '\n') {
                if (word.length() > 0) {
                    word += (char) 0;
                    ++lib[word];
                    word = "";
                }
            } else {
                word += s[i];
            }
        }
    }
    return lib;
}

int main(int argc, const char * argv[]) {
    const int n = 13;
    std::cout << "Fibonaci number " << n << " = " << Fibonacci<n>::f << std::endl;
    std::map<std::string, int> lib = mapper("/Users/Viktor/Documents/MIPT_base/prog/homework_december/homework_december/file.txt")  ;
    std::map<std::string, int>::iterator iter;
    for (iter = lib.begin(); iter != lib.end(); iter++)
        std::cout << iter->second << "\t| " << iter->first << std::endl;
    return 0;
}
