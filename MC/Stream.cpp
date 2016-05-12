#include "Stream.h"

namespace MC {

Stream::Stream() { reset(); }

void Stream::reset() {
  showHeader(true);
  dAction(true);
  dStateStep(1);
  dStatStep(1);
  dInfoStep(1);
  actionFile("action.log");
  stateFile("state.log");
  statFile("stat.log");
  infoFile("info.log");
  faction.clear();
  fstate.clear();
  fstat.clear();
  finfo.clear();
}

void Stream::setWorkDir(const std::string& _wdir) {
  actionFile() = _wdir+"/"+actionFile();
  stateFile() = _wdir+"/"+stateFile();
  statFile() = _wdir+"/"+statFile();
  infoFile() = _wdir+"/"+infoFile();
}

bool Stream::init(int _precision, std::ios::openmode _mode, bool _silent) {
  faction.open(actionFile().c_str(),_mode);
  fstate.open(stateFile().c_str(),_mode);
  fstat.open(statFile().c_str(),_mode);
  finfo.open(infoFile().c_str(),_mode);

  setFormat(_precision);

  return check(_silent);
}

void Stream::final() {
  faction.close();
  fstate.close();
  fstat.close();
  finfo.close();
}

std::ofstream& Stream::action(bool auto_width) {
  if(auto_width) faction.width(int(faction.precision())+9);
  return faction;
}

std::ofstream& Stream::state(bool auto_width) {
  if(auto_width) fstate.width(int(fstate.precision())+9);
  return fstate;
}

std::ofstream& Stream::stat(bool auto_width) {
  if(auto_width) fstat.width(int(fstat.precision())+9);
  return fstat;
}

std::ofstream& Stream::info(bool auto_width) {
  if(auto_width) finfo.width(int(finfo.precision())+9);
  return finfo;
}

void Stream::setFormat(int _precision) {
  faction.setf(std::ios::scientific, std::ios::floatfield);
  fstate.setf(std::ios::scientific, std::ios::floatfield);
  fstat.setf(std::ios::scientific, std::ios::floatfield);
  finfo.setf(std::ios::scientific, std::ios::floatfield);
  faction.precision(_precision);
  fstate.precision(_precision);
  fstat.precision(_precision);
  finfo.precision(_precision);
  faction << std::left;
  fstate << std::left;
  fstat << std::left;
  finfo << std::left;
}

bool Stream::check(bool _silent) {
  if(!_silent) {
    if(!faction.good()) {
      std::cout << "WARNING: Cannot open `" 
		<< actionFile() << "'" << std::endl; 
    }
    if(!fstate.good()) {
      std::cout << "WARNING: Cannot open `" 
		<< stateFile() << "'" << std::endl; 
    }
    if(!fstat.good()) {
      std::cout << "WARNING: Cannot open `" 
		<< statFile() << "'" << std::endl; 
    }
    if(!finfo.good()) {
      std::cout << "WARNING: Cannot open `" 
		<< infoFile() << "'" << std::endl; 
    }
  }
  return (faction.good()&&fstate.good()&&fstat.good()&&finfo.good());
}

}/* MC */
