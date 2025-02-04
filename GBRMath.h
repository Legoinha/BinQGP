#ifndef GBRMATH
#define GBRMATH

#include "RVersion.h"
#if ROOT_VERSION_CODE >=  ROOT_VERSION(5,34,00)
  #include "vdtMath.h"
#else
  #include "DataFormats/Math/interface/VDTMath.h" 
#endif
  
#include <cmath>

namespace gbrmath {
 
  inline double fast_pow(double base, double exponent) {
    if (base==0. && exponent>0.) return 0.;
    else if (base>0.) return vdt::fast_exp(exponent*vdt::fast_log(base));
    else return std::nan("");
  }

//   inline float fast_powf(float base, float exponent) {
//     if (base==0. && exponent>0.) return 0.;
//     else if (base>0.) return vdt::fast_expf(exponent*vdt::fast_logf(base));
//     else return std::nanf(""); 
//   }

}

#endif
