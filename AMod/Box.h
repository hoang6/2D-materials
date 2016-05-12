#ifndef _AMOD_BOX_H
#define _AMOD_BOX_H

#include <vector>

namespace AMod {

class Box {
public:
  struct Information {
    int chiV[2];
    int heightV[2];
    int iQuad[2][2];
    double dQuad[2][2];
    double chi;
    int colRange[2];
    int rowRange[2];
    int nGrids;
  };
  //neighbor positions
  static const int TOP_LEFT;
  static const int TOP_RIGHT;
  static const int BOTTOM_LEFT;
  static const int BOTTOM_RIGHT;
  //scale
  static const double SCALE;

  Box();
  Box(const Box& box);
  Box& operator= (const Box& box);
  void prepare(const int _chiV[2], const int _heightV[2]);
  bool inQuad(int c, int r) const;
  void move2quad(int c, int r, int& cp, int& rp) const;
  static void neighbor(int c, int r, int pos, int& neighb_c, int& neighb_r);
  static void cr2xy(int c, int r, double& x, double& y);
  static void xy2cr(double x, double y, int& c, int& r);
  static double length(int c, int r);
  const Box::Information& information() const;

private:
  Box::Information info;
  std::vector<int> ranges[2];

  bool onBoundary(int c, int r) const;
  double computeRange(int r, int id) const;
  void computeRanges();
  static double distance2(double x1, double y1, double x2, double y2);
};

/************************************************************/
inline const Box::Information& Box::information() const { return info; }

} /* AMOD */

#endif
