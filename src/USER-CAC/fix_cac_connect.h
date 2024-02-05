#ifdef FIX_CLASS

FixStyle(cac/connect, FixConnectCAC)

#else

#ifndef LMP_FIX_CONNECT_CAC_H
#define LMP_FIX_CONNECT_CAC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixConnectCAC : public Fix {
 public:
  FixConnectCAC(class LAMMPS *, int, char **);
  virtual ~FixConnectCAC() {}
  int setmask();
  virtual void setup(int);

 protected:
};

}    // namespace LAMMPS_NS

#endif
#endif
