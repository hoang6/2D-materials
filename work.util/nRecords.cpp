#include <string>
#include <fstream>
#include <iostream>

int main(int argc, char* argv[]) {
  if(argc != 3) {
    std::cout << "Usage: nRecords <file> <^key>" << std::endl;
    return 1;
  }

  std::string fname = argv[1];
  std::string key = argv[2];

  std::fstream fin(fname.c_str());
  if(!fin.good()) {
    std::cout << "Error: cannot open file: " << fname << std::endl;
    return 2;
  }
  std::string tline;
  int count = 0;
  while(!fin.eof()) {
    std::getline(fin,tline);
    if(tline.substr(0,key.length()) == key) count++;
  }
  fin.close();

  std::cout << "Read " << count << " records" << std::endl;

  return 0;
}
