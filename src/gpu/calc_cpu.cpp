//compute_cpu.cpp

#include <vector>
#include "CALC_CPU_H_cpu.h"
template<typename T>

computeCellAverageTemperature()


computeLocalCapacityMatrix()
{/*
 T avgTemperature = localAverageTemperature[];//precompute avrgae temperature with shape functions
 T avgFactor = localDensity[] * getCapacity(avgTemperature) * localWeight[] * std::abs(localJacobian[]);
  for(int i = 0; i < supportNodeSize; ++i){
    for(int j = 0; j < supportNodeSize; ++j){
      C.writeElement( avgFactor * LocalShapeFunPhi[(0,i)] * LocalShapeFunPhi[(0,i)], i, j );
    }
  }*/
}

computeLocalConductivityMatrix()
{/*
  int max_deriv_index = dim + 1;
  T avgTemperature = localAverageTemperature[];//precompute average temperature with shape functions
  T avgFactor = mat->getKappa(avgTemp) * weight * std::abs(jacobian);
   for(int i = 0; i < supportNodeSize; ++i){
     for(int j = 0; j < supportNodeSize; ++j){
       for (auto m = 1; m < max_deriv_index; ++m) {
             H.addElement( (LocalShapeFunPhi[(0,i)] * LocalShapeFunPhi[(0,i)] * avgFactor,i,j);
       }
     }
   }
*/
}

/*
double interpolate1D(double key, const std::map<double, double>& theMap)
{
    typedef std::map<double, double>::const_iterator i_t;

    i_t i = theMap.upper_bound(key);
    if (i == theMap.end()) {
        return (--i)->second;
    }
    if (i == theMap.begin()) {
        return i->second;
    }
    i_t l = i;
    --l;

    const double delta = (key - l->first) / (i->first - l->first);
    return (delta * i->second + (1 - delta) * l->second);
}
*/
