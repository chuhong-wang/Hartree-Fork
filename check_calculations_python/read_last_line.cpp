#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
// #include <filesystem>

using std::string, std::ifstream, std::cerr, std::istringstream, std::cout, std::endl;

auto main() -> int {
  int num = 5;
  string filename = "../data/S.DAT";

  // cout << "Full Path: " << filesystem::absolute(filename) << endl;

  int nao;
  string firstline;
  ifstream fin(filename);
  if (fin.fail()) {
    cerr << "file doesn't exist \n";
  }

  std::string line;
  std::string last_line;

  // Read the file line by line
  while (std::getline(fin, line)) {
    last_line = line;
  }

  // Close the file
  fin.close();

  istringstream sin(last_line);
  sin >> nao;
  cout << nao << '\n';

  return 0;
}