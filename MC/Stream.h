#ifndef _MC_STREAM_H
#define _MC_STREAM_H

#include "../Util/Accessors.h"
#include <string>
#include <fstream>
#include <iostream>

namespace MC {

class Stream {
public:
  Stream();
  void reset();
  void setWorkDir(const std::string& _wdir);
  bool init(int _precision, std::ios::openmode _mode, bool _silent);
  void final();
  std::ofstream& action(bool auto_width = true);
  std::ofstream& state(bool auto_width = true);
  std::ofstream& stat(bool auto_width = true);
  std::ofstream& info(bool auto_width = true);
  Util::Accessors<bool> showHeader;
  Util::Accessors<bool> dAction;
  Util::Accessors<long> dStateStep;
  Util::Accessors<long> dStatStep;
  Util::Accessors<long> dInfoStep;
  Util::Accessors<std::string> actionFile;
  Util::Accessors<std::string> stateFile;
  Util::Accessors<std::string> statFile;
  Util::Accessors<std::string> infoFile;

private:
  void setFormat(int _precision);
  bool check(bool _silent);
  Stream(const Stream& stream);
  Stream& operator= (const Stream& stream);

  std::ofstream faction;
  std::ofstream fstate;
  std::ofstream fstat;
  std::ofstream finfo;
};

}/* MC */

#endif
