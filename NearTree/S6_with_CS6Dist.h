#ifndef S6_WITH_CS6DIST_H
#define S6_WITH_CS6DIST_H
#include "S6.h"
/* #include "CS6Dist.h" */
double CS6Dist(double ds1[6], double ds2[6]);
#include <cmath>
class S6_with_CS6Dist : public S6
{
  public: 
    using S6::S6;
    double DistanceBetween( const S6& s2 ) {
      double ds1[6], ds2[6];
      size_t ii;
      for (ii=0; ii < 6; ii++) {
        ds1[ii]=(*this)[ii]; ds2[ii]=s2[ii];
      }
      return 0.1*std::sqrt(CS6Dist(ds1,ds2));
   }
};
#endif

