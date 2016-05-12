#ifndef _AMOD_MOL_ADJUSTER_H
#define _AMOD_MOL_ADJUSTER_H

#include "MolService.h"
#include "../Util/Array1D.h"
#include "../Util/Array2D.h"
#include <utility>
#include <vector>

namespace AMod {

class MolAdjuster: public MolServ {
public:
  typedef Util::Array1D<double,3> Point;
  typedef Util::Array2D<double,3,3> Axes;

  MolAdjuster();
  MolAdjuster(Molecule& _mol);
  void reset();
  /****** mass center ******/
  void adjustMassCenter(const Point& _massCenter = Point(0.0)) const;
  Point massCenter() const;
  static Point massCenter(const Molecule& _mol);
  /****** principle axes ******/
  void adjustPrinpAxes(const Axes& _prinpAxes = Axes(1.0,0.0)) const;
  Axes prinpAxes() const;
  static Axes prinpAxes(const Molecule& _mol);
  /****** adjust mass center and prinp axes ******/
  void adjust(const Point& _massCenter = Point(0.0), 
	      const Axes& _prinpAxes = Axes(1.0,0.0)) const;
  void adjust(const Molecule& _mol) const;
  /****** move & rot ******/
  void move(const Point& _oldPt, const Point& _newPt) const;
  void rotate(const Axes& _oldAxes, const Axes& _newAxes) const;
  /****** align ******/
  static bool align(Molecule& mol, int bondID, 
		    const Molecule& _mol, int _bondID);

protected:
  static void inertial(const Molecule& _mol,
		       double& Ixx, double& Iyy, double& Izz,
		       double& Ixy, double& Iyz, double& Izx);
  static void adjustPosition(Molecule& _mol,
			     const Point& _oldPt, const Point& _newPt);
  static void adjustAxes(Molecule& _mol, 
			 const Axes& _oldAxes, const Axes& _newAxes);
};

/************************************************************/
inline MolAdjuster::Point MolAdjuster::massCenter() const { return massCenter(molecule()); }

inline MolAdjuster::Axes MolAdjuster::prinpAxes() const { return prinpAxes(molecule()); }

} /* AMod */

#endif
