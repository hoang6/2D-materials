#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 2) {
    cout << "Usage: simpleMPIDriver <cmd_file>" << endl;
    return 1;
  }

  ifstream fin(argv[1]);
  if(!fin.good()) {
    cout << "ERROR: Cannot open file " << argv[1] << endl;
    return 2;
  }
  string cmd;
  vector<string> cmds;
  while(getline(fin,cmd)) {
    if(!cmd.empty())
      cmds.push_back(cmd);
  }
  fin.close();

  MPI::Init();

  int ntasks = cmds.size();
  int np = MPI::COMM_WORLD.Get_size();
  int id = MPI::COMM_WORLD.Get_rank();
  int chunk = int(ceil(double(ntasks)/np));
  int start = chunk*id;
  int end = (id == np-1 ? ntasks-1 : chunk*id+chunk-1);
  
  for(int k = start; k <= end; k++) {
    cout << "This is id " << id << " doing task " << k << endl;
    cout << cmds[k] << endl;
    system(cmds[k].c_str());
  }

  MPI::Finalize();
  
  return 0;
}
