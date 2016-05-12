#include "../AnyOption/ParseMainArgs.h"
#include "../AMod/AMod.h"
#include "../Util/util.h"
#include <cmath>

int fixCenter(AMod::Molecule& _mol);
void fixCNTPBC(AMod::Molecule& _mol);

int main(int argc, char* argv[]) {
  /********** Parse arguments **********/
  ParseMainArgs::Input input;
  std::ostringstream ossUsage;
  ParseMainArgs* par = new ParseMainArgs;
  (*par) << "Usage: manufacture [OPTION]... SIZE"
	 << "Manufacture a carbon molecule"
	 << ""
	 << " --help                -h       print this help"
	 << " --mol                 -m [s]   type of molecule"
	 << "                                C60 : fullerene"
	 << "                                GPH : graphene"
	 << "                                CNT : carbon nanotube"
	 << " --NPBC                -N       none periodic boundary condition"
	 << " --XH                  -X       no capping hydrogen atoms"
	 << " --fix_gesture         -g       fix the gesture of the molecule"
	 << " --auto_file           -a       output to an auto named file"
	 << " <'i i i i'>                    chirality (1,2), height(1,2)"
  	 << "                                hint: |_| => m'/n'=-(m+2n)/(2m+n)";
  par->process("help h|NPBC N|XH X|fix_gesture g|auto_file a", "mol m", argc, argv, input);
  par->printUsage(ossUsage);
  delete par;

  /********** Read arguemnts **********/
  if(input.flags.find("help h")->second) {
    std::cout << ossUsage.str();
    return 1;
  }

  int chiV[2];
  int heightV[2];
  std::string molType = input.options.find("mol m")->second;
  bool fixGesture = input.flags.find("fix_gesture g")->second;

  if(molType == "GPH" || molType == "CNT") {
    if(input.args.size() != 1) {
      std::cout << "manufacture: molecule size is not spectified" << std::endl;
      return 2;
    }
    ParseMainArgs::StringArray words;
    ParseMainArgs::split(input.args[0], ' ', words);
    if(words.size() != 4 ||
       !Util::readWord(words[0], chiV[0]) || 
       !Util::readWord(words[1], chiV[1]) ||
       (chiV[0] == 0 && chiV[1] == 0) ||
       !Util::readWord(words[2], heightV[0]) ||
       !Util::readWord(words[3], heightV[1]) ||
       (heightV[0] == 0 && heightV[0] == 0)) {
      std::cout << "manufacture: illegal molecule size" << std::endl;
      return 3;
    }
  }

  /********** Manufacture **********/
  AMod::Molecule mol;
  AMod::MolFactory mfactory(mol);
  std::ostringstream ossFile;

  if(molType == "C60") {//make C60
    mfactory.makeC60();
    ossFile << "C60_";
    if(fixGesture) ;
  }
  else if(molType == "GPH") {
    ossFile << "GPH_";
    if(!input.flags.find("NPBC N")->second) {//make GPH_PBC
      mfactory.makeGPH_PBC(chiV, heightV);
      ossFile << "PBC_";
      if(fixGesture) fixCenter(mol);
    }
    else {//make GPH_NPBC
      mfactory.makeGPH_NPBC(chiV, heightV, !input.flags.find("XH X")->second);
      ossFile << "NPBC_";
      if(fixGesture) ;
    }
  }
  else if(molType == "CNT") {
    ossFile << "CNT_";
    if(!input.flags.find("NPBC N")->second) {//make CNT_PBC
      mfactory.makeCNT_PBC(chiV, heightV);
      ossFile << "PBC_";
      if(fixGesture) fixCNTPBC(mol);
    }
    else {//make CNT_NPBC
      mfactory.makeCNT_NPBC(chiV, heightV, !input.flags.find("XH X")->second);
      ossFile << "NPBC_";
      if(fixGesture) ;
    }
  }
  else {
    std::cout << "manufacture: unrecognized mol type: " << molType << std::endl;
    std::cout << "Try `manufacture --help' for more information" << std::endl;
    return 4;
  }

  for(int i = 0; i < int(mol.naxes()); i++)
    for(int j = 0; j < 3; j++)
      if(std::abs(mol.axis(i).x[j]) < Util::EPS_DOUBLE)
	mol.axis(i).x[j] = 0.0;

  if(mol.natoms() == 0) {
    std::cout << "manufacture: empty molecule" << std::endl;
    return 5;
  }

  ossFile << "C" << mol.countAtoms(AMod::PTE::Carbon) 
	  << "H" << mol.countAtoms(AMod::PTE::Hydrogen);

  if(input.flags.find("auto_file a")->second) {
    std::ofstream fout(ossFile.str().c_str());
    AMod::MolO::setFormat(fout);
    if(!fout.good()) {
      std::cout << "manufacture: cannot open file: " << ossFile.str() << std::endl;
      return 6;
    }
    std::cout << "manufacture: produce file: " << ossFile.str() << std::endl;
    mol.io().dumpTxt(fout);
    fout.close();
  }
  else {
    AMod::MolO::setFormat(std::cout);
    mol.io().dumpTxt(std::cout);
  }

  return 0;
}

int fixCenter(AMod::Molecule& _mol) {
  AMod::MolAdjuster madj(_mol);
  int atomID = _mol.topo().atomID(madj.massCenter().data);
  for(int k = 0; k < 3; k++) _mol.atom(atomID)->fixed[k] = 1;
  return atomID;
}

void fixCNTPBC(AMod::Molecule& _mol) {
  int atomID = fixCenter(_mol);
  const AMod::Atom& atom = _mol.atom(atomID);

  int k, kmax = atom.arrows.size();
  int kNeighb = -1;
  double score = 0.0, tscore = 0.0;
  double dx[3], tdx[3];
  for(k = 0; k < kmax; k++) {
    Util::vcopy(atom.arrows[k].dx,tdx,3);
    Util::normalize(tdx,3);
    tscore = Util::abs(tdx[2]);
    if(kNeighb < 0 || score > tscore) {
      kNeighb = k;
      score = tscore;
      Util::vcopy(tdx,dx,3);
    }
  }
  if(kNeighb < 0) return;
  
  AMod::MolAdjuster madj(_mol);
  double th = Util::PI/2.0-std::atan2(dx[1],dx[0]);
  AMod::MolAdjuster::Axes oldAxes(1.0,0.0), newAxes(1.0,0.0);
  newAxes[0][0] = +std::cos(th); newAxes[0][1] = +std::sin(th);
  newAxes[1][0] = -std::sin(th); newAxes[1][1] = +std::cos(th);
  madj.rotate(oldAxes,newAxes);
  atom.arrows[kNeighb].moon()->fixed[0] = 1;
}
