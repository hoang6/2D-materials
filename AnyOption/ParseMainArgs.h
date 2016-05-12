#ifndef _PARSE_MAIN_ARGS_H
#define _PARSE_MAIN_ARGS_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

class AnyOption;
class ParseMainArgs {
public:
  typedef std::vector<bool> BoolArray;
  typedef std::vector<std::string> StringArray;
  class Input;

  ParseMainArgs();
  ~ParseMainArgs();
  ParseMainArgs& operator<< (const char* line);
  void printUsage(std::ostream& _fout = std::cout);
  void process(const std::string& flag_names,   //e.g. "help h|flag1 f"
	       const std::string& option_names,
	       int argc, char* argv[], Input& input);
  static void echo(Input& input);
  static int nSegments(const std::string& line);
  static void split(const std::string& line, StringArray& segments);
  static int nSegments(const std::string& line, char delim);
  static void split(const std::string& line, char delim, StringArray& segments);

  struct StrComp {
    bool operator() (const std::string& lhs, const std::string& rhs) const {
      return lhs.compare(rhs) < 0; } 
  };
  typedef std::map<std::string,bool,StrComp> Flags;
  typedef std::map<std::string,std::string,StrComp> Options;
  class Input { 
  public:
    Flags flags; 
    Options options; 
    StringArray args; 
    
    Input() { reset(); }
    void reset() { flags.clear(); options.clear(); args.clear(); }
  };

private:
  static void split2pairs(const std::string& line, StringArray& segments);

  AnyOption* opt;
};

#endif
