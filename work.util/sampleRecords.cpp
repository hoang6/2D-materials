#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "../Util/util.h"

int main(int argc, char* argv[]) {
  if(argc != 6) {
    std::cout << "Usage: sampleRecords <file> <^key> <start> <interval> <N>" 
	      << std::endl;
    return 1;
  }

  std::string fname = argv[1];
  std::string key = argv[2];
  int a, d, n;
  if(!Util::readWord(argv[3],a) || a == 0) {
    std::cout << "Error: illegal <start>" << std::endl;
    return 2;
  }
  if(!Util::readWord(argv[4],d) || d <= 0) {
    std::cout << "Error: illegal <interval>" << std::endl;
    return 3;
  }
  if(!Util::readWord(argv[5],n) || n <= 0) {
    std::cout << "Error: illegal <N>" << std::endl;
    return 4;
  }

  std::fstream fin(fname.c_str());
  if(!fin.good()) {
    std::cout << "Error: cannot open file: " << fname << std::endl;
    return 5;
  }
  std::string tline;
  int count;
  std::ifstream::pos_type posStart;
  if(a > 0) {
    count = 0;
    while(!fin.eof()) {
      posStart = fin.tellg();
      std::getline(fin,tline);
      if(tline.substr(0,key.length()) == key) count++;
      if(count == a) { fin.seekg(posStart); fin.clear(); break; }
    }
  }
  else {
    char buf;
    std::ifstream::pos_type posBeg = fin.tellg();
    fin.seekg(-1,std::ios::end);
    count = 0;
    while(true) {
      buf = static_cast<char>(fin.peek());
      if(fin.tellg() == posBeg) tline += buf;
      if(buf == '\n' || fin.tellg() == posBeg) {
	std::reverse(tline.begin(),tline.end());
	if(tline.substr(0,key.length()) == key) count++;
	if(count == -a) {
	  if(fin.tellg() != posBeg) fin.seekg(+1,std::ios::cur);
	  fin.clear();
	  break;
	}
	tline.clear();
      }
      else tline += buf;
      if(fin.tellg() == posBeg) {
	fin.seekg(0,std::ios::end);
	break;
      }
      fin.seekg(-1,std::ios::cur);
    }
  }
  tline.clear();
  count = 0;
  while(!fin.eof()) {
    std::getline(fin,tline);
    if(tline.substr(0,key.length()) == key) count++;
    if(count > n) { count--; break; }
    if(count > 0 && (count-1)%d == 0)
      std::cout << tline << std::endl;
  }
  fin.close();

  std::cerr << "Read " << (count-1)/d+1 << " records" << std::endl;

  return 0;
}
