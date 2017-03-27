//some helper functions, arrangeds to fulfill Structure Of Arrays (SOA) philoshophy
#include "functions_cpu.h"

double computeCellAverageTemperature(){
  double avg=0.0;
  return avg;
}

double getPhi_from_table(ShapeFunctionTable &shapeFunctionTable,
                         int gp_id,
                         int derivative_order,
                         int node_id)
{
  int ndim = shapeFunctionTable.dim;
  int num_nodes = shapeFunctionTable.number_points;
  return shapeFunctionTable.phis[node_id * (ndim * num_nodes) + derivative_order];
}

void setPhi_into_table(double phiValue,
                       ShapeFunctionTable &shapeFunctionTable,
                       int gp_id,
                       int derivative_order,
                       int node_id)
{
  int ndim = shapeFunctionTable.dim;
  int num_nodes = shapeFunctionTable.number_points;
  shapeFunctionTable.phis[node_id * (ndim * num_nodes) + derivative_order] = phiValue;
}

/*double interpolate1D(double key, const std::map<double, double>& theMap)
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
}*/
double interpolate1D(double query_value,
                    double *reference_values,
                    double *sought_values,
                    int init_position,
                    int counter)
{
  bool upper_bounded = false;
  int upper_index = 0;
  for(int i = 0; i < counter; i++){
    //first find the upper_bound
    double this_val = reference_values[init_position + counter];
    if(query_value <= this_val){//bounded here
      upper_index = i;
      if(upper_index == 0) return sought_values[init_position];
      upper_bounded = true;
    }
  }
  if(!upper_bounded) return sought_values[init_position + counter];//not bound found so return last
  //so we have a bound
  int lower_index = upper_index - 1;
  double delta = (query_value - reference_values[init_position + lower_index]) / (reference_values[init_position + upper_index] - reference_values[init_position + lower_index]);
  return delta * sought_values[init_position + upper_index] + (1.0 - delta) * sought_values[init_position + lower_index];
}


void addMaterialToTable(){}

double getMaterialKappa (MaterialTable &materials,
                    int material_id,
                    double average_temperature)
{
  int n_vals =  materials.m_kappa_counters[material_id];
  int init_vals =  materials.m_kappa_inits[material_id];
  return interpolate1D(average_temperature,
                       materials.m_kappa_temps,
                       materials.m_kappa_values,
                       init_vals,
                       n_vals);
}

double getMaterialDensity (MaterialTable &materials,
                      int material_id)
{
  return materials.density[material_id];
}

double getMaterialCapacity (MaterialTable &materials,
                      int material_id,
                      double average_temperature)
{
  int n_vals =  materials.m_capacity_counters[material_id];
  int init_vals =  materials.m_capacity_inits[material_id];
  return interpolate1D(average_temperature,
                       materials.m_capacity_temps,
                       materials.m_capacity_values,
                       init_vals,
                       n_vals);

}
