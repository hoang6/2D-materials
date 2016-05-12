#include "AnyOption.h"
#include "ParseMainArgs.h"

ParseMainArgs::ParseMainArgs() {
  opt = new AnyOption;
}

ParseMainArgs::~ParseMainArgs() {
  delete opt;
}

ParseMainArgs& ParseMainArgs::operator<< (const char* line) {
  opt->addUsage(line);
  return *this;
}

void ParseMainArgs::printUsage(std::ostream& _fout) { opt->printUsage(_fout); }

void ParseMainArgs::process(const std::string& flag_names,
			    const std::string& option_names,
			    int argc, char* argv[], Input& input) {
  StringArray flag_na;
  StringArray option_na;
  bool valid[2];
  std::string key;
  Flags& flags = input.flags;
  Options& options = input.options;
  StringArray& args = input.args;

  flags.clear();
  options.clear();
  args.clear();

  /*** preferences ***/
  //opt->noPOSIX(); /* do not check for POSIX style character options */
  //opt->setVerbose(); // print warnings about unknown options
  //opt->autoUsagePrint(true); // print usage for bad options
  /*** set flags ***/
  split2pairs(flag_names, flag_na);
  for(int k = 0, kmax = flag_na.size(); k < kmax; k += 2) {
    valid[0] = bool(flag_na[k  ].size());
    valid[1] = bool(flag_na[k+1].size());
    if(valid[0] && valid[1])
      opt->setCommandFlag(flag_na[k].c_str(),flag_na[k+1][0]);
    else if(valid[0] && !valid[1])
      opt->setCommandFlag(flag_na[k].c_str());
    else if(!valid[0] && valid[1])
      opt->setCommandFlag(flag_na[k+1][0]);
  }
  /*** set options ***/
  split2pairs(option_names, option_na);
  for(int k = 0, kmax = option_na.size(); k < kmax; k += 2) {
    valid[0] = bool(option_na[k  ].size());
    valid[1] = bool(option_na[k+1].size());
    if(valid[0] && valid[1])
      opt->setCommandOption(option_na[k].c_str(),option_na[k+1][0]);
    else if(valid[0] && !valid[1])
      opt->setCommandOption(option_na[k].c_str());
    else if(!valid[0] && valid[1])
      opt->setCommandOption(option_na[k+1][0]);
  }
  /*** process argmuments ***/
  opt->processCommandArgs(argc, argv);
  /*** get flags ***/
  for(int k = 0, kmax = flag_na.size(); k < kmax; k += 2) {
    valid[0] = bool(flag_na[k  ].size());
    valid[1] = bool(flag_na[k+1].size());
    if(valid[0] && valid[1]) {
      key = flag_na[k]+" "+flag_na[k+1];
      flags[key] = (opt->getFlag(flag_na[k].c_str()) || 
		    opt->getFlag(flag_na[k+1][0]));
    }
    else if(valid[0] && !valid[1]) {
      key = flag_na[k];
      flags[key] = opt->getFlag(flag_na[k].c_str());
    }
    else if(!valid[0] && valid[1]) {
      key = flag_na[k+1];
      flags[key] = opt->getFlag(flag_na[k+1][0]);
    }
  }
  /*** get options ***/
  for(int k = 0, kmax = option_na.size(); k < kmax; k += 2) {
    valid[0] = bool(option_na[k  ].size());
    valid[1] = bool(option_na[k+1].size());
    if(valid[0] && valid[1]) {
      key = option_na[k]+" "+option_na[k+1];
      if(opt->getValue(option_na[k].c_str()) != NULL)
	options[key] = opt->getValue(option_na[k].c_str());
      else if(opt->getValue(option_na[k+1][0]) != NULL)
	options[key] = opt->getValue(option_na[k+1][0]);
      else
	options[key] = "";
    }
    else if(valid[0] && !valid[1]) {
      key = option_na[k];
      if(opt->getValue(option_na[k].c_str()) != NULL)
	options[key] = opt->getValue(option_na[k].c_str());
      else
	options[key] = "";
    }
    else if(!valid[0] && valid[1]) {
      key = option_na[k+1];
      if(opt->getValue(option_na[k+1][0]) != NULL)
	options[key] = opt->getValue(option_na[k+1][0]);
      else
	options[key] = "";
    }
  }
  /*** other arguments ***/
  for(int k = 0; k < opt->getArgc(); k++) 
    args.push_back(std::string(opt->getArgv(k)));
}

int ParseMainArgs::nSegments(const std::string& line) {
  StringArray segments;
  split(line, segments);
  return segments.size();
}

void ParseMainArgs::split(const std::string& line, StringArray& segments) {
  std::istringstream iss(line);
  std::string s;
  while(iss >> s) segments.push_back(s);
}

int ParseMainArgs::nSegments(const std::string& line, char delim) {
  StringArray segments;
  split(line, delim, segments);
  return segments.size();
}

void ParseMainArgs::split(const std::string& line, char delim, StringArray& segments) {
  std::istringstream iss(line);
  std::string s;
  while(std::getline(iss,s,delim)) segments.push_back(s);
}

void ParseMainArgs::split2pairs(const std::string& line, StringArray& segments) {
  StringArray sa, sub_sa;
  split(line, '|', sa);
  for(int k = 0, kmax = sa.size(); k < kmax; k++) {
    split(sa[k], ' ', sub_sa);
    if(sub_sa.size() == 1) {
      if(sub_sa[0].size() > 1) {//...|flag|...
	segments.push_back(sub_sa[0]);
	segments.push_back("");
      }
      else if(sub_sa[0].size() == 1) {//...|f|...
	segments.push_back("");
	segments.push_back(sub_sa[0]);
      }
    }
    else if(sub_sa.size() > 1) {
      //...|... flag ...|...
      segments.push_back("");
      for(int j = 0, jmax = sub_sa.size(); j < jmax; j++) {
	if(sub_sa[j].size() > 1) {
	  segments.back() = sub_sa[j];
	  break;
	}
      }
      //...|... f ...|...
      segments.push_back("");
      for(int j = 0, jmax = sub_sa.size(); j < jmax; j++) {
	if(sub_sa[j].size() == 1) {
	  segments.back() = sub_sa[j];
	  break;
	}
      }
    }
    sub_sa.clear();
  }//for(int k...
}
