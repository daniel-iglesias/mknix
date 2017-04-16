//some helper functions, arrangeds to fulfill Structure Of Arrays (SOA) philoshophy
#include "functions_cpu.h"
#include <stdlib.h>
#include <map>

double getPhi_from_table(ShapeFunctionTable *shapeFunctionTable,
                         int point_id,
                         int local_node_id)
{
  int support_node_size = shapeFunctionTable->support_node_size;
  return shapeFunctionTable->phis[point_id * support_node_size + local_node_id];
}

void setPhi_into_table(double phiValue,
                       ShapeFunctionTable *shapeFunctionTable,
                       int point_id,
                       int local_node_id)
{
  int support_node_size = shapeFunctionTable->support_node_size;
  shapeFunctionTable->phis[point_id * support_node_size + local_node_id] = phiValue;
}

void init_shape_functions_table(ShapeFunctionTable **shapeFunctionTable,
                                int support_node_size,
                                int number_points)
{
  (*shapeFunctionTable)->support_node_size = support_node_size;
  (*shapeFunctionTable)->number_points = number_points;
  (*shapeFunctionTable)->phis = (double*) malloc(number_points * support_node_size *  sizeof(double));
}

void free_shape_functions_table(ShapeFunctionTable **shapeFunctionTable)
{
  free((*shapeFunctionTable)->phis);
}

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

double getMaterialKappa (MaterialTable *materials,
                        int material_id,
                        double average_temperature)
{
  int n_vals =  materials->_kappa_counters[material_id];
  int init_vals =  materials->_kappa_inits[material_id];
  return interpolate1D(average_temperature,
                       materials->_kappa_temps,
                       materials->_kappa_values,
                       init_vals,
                       n_vals);
}

double getMaterialDensity (MaterialTable *materials,
                      int material_id)
{
  return materials->density[material_id];
}

double getMaterialCapacity (MaterialTable *materials,
                            int material_id,
                            double average_temperature)
{
  int n_vals =  materials->_capacity_counters[material_id];
  int init_vals =  materials->_capacity_inits[material_id];
  return interpolate1D(average_temperature,
                       materials->_capacity_temps,
                       materials->_capacity_values,
                       init_vals,
                       n_vals);

}

void setupMaterialTables(MaterialTable **mt,
                         std::map<int, mknix::Material> &materials)
{
    //MaterialTable *mt;
    /*for (auto& system : subSystems) {
        mt = system.second->getMaterialTablePointer();
    }*/
    std::cout << "setting up SOA material tables" << std::endl;
    int nMat = materials.size();
    std::cout << "Allocating memory for "<< nMat << " materials" << std::endl;
    (*mt)->number_materials = nMat;
    (*mt)->capacity = (double*) malloc(nMat * sizeof(double));
    (*mt)->kappa = (double*) malloc(nMat * sizeof(double));
    (*mt)->beta = (double*) malloc(nMat * sizeof(double));
    (*mt)->density = (double*) malloc(nMat * sizeof(double));
    (*mt)->_capacity_counters = (int*) malloc(nMat * sizeof(int));
    (*mt)->_capacity_inits = (int*) malloc(nMat * sizeof(int));
    (*mt)->_kappa_counters = (int*) malloc(nMat * sizeof(int));
    (*mt)->_kappa_inits = (int*) malloc(nMat * sizeof(int));
    int i = 0;
    std::cout << "Copy variables for materials" << std::endl;
    for(auto const &imat : materials){
      (*mt)->density[i] = imat.second.retDensity();
      (*mt)->beta[i] = imat.second.retBeta();
      (*mt)->capacity[i] = imat.second.retCapacity();
      (*mt)->_capacity_counters[i]= imat.second.retCapacityMapSize();
      (*mt)->kappa[i] = imat.second.retKappa();
      (*mt)->_kappa_counters[i]= imat.second.retKappaMapSize();
      i++;
    }
    i = 0;
    for(auto const &imat : materials){
      (*mt)->capacity[i] = imat.second.retCapacity();
      (*mt)->_capacity_counters[i]= imat.second.retCapacityMapSize();
      (*mt)->kappa[i] = imat.second.retKappa();
      (*mt)->_kappa_counters[i]= imat.second.retKappaMapSize();
      i++;
    }
    (*mt)->_capacity_inits[0] = 0;
    (*mt)->_kappa_inits[0] = 0;
    int capacity_tot_vals = (*mt)->_capacity_counters[0];
    int kappa_tot_vals = (*mt)->_kappa_counters[0];
    for(int imat = 1; imat < nMat; imat++){//prefix sum scan
      capacity_tot_vals = (*mt)->_capacity_counters[imat];
      kappa_tot_vals = (*mt)->_kappa_counters[imat];
      (*mt)->_capacity_inits[imat] = (*mt)->_capacity_inits[imat - 1] + (*mt)->_capacity_counters[imat-1];
      (*mt)->_kappa_inits[imat] = (*mt)->_kappa_inits[imat - 1] + (*mt)->_kappa_counters[imat-1];
    }
    (*mt)->_capacity_temps = (double*) malloc(capacity_tot_vals * sizeof(double));
    (*mt)->_capacity_values = (double*) malloc(capacity_tot_vals * sizeof(double));
    (*mt)->_kappa_temps = (double*) malloc(kappa_tot_vals * sizeof(double));
    (*mt)->_kappa_values = (double*) malloc(kappa_tot_vals* sizeof(double));
    //now we copy all values
    i = 0;
    for(auto const &imat : materials){
      const std::map<double, double> *cap_map = imat.second.retCapacityMap();
      int initc = (*mt)->_capacity_inits[i];
      //for(int icap = 0;  icap < cap_map->size(); icap++){
      int j =0;
      //for (auto const& icap:cap_map){
      for(std::map<double,double>::const_iterator icap = cap_map->begin(); icap != cap_map->end(); ++icap){
        (*mt)->_capacity_temps[initc + j] = icap->first;
        (*mt)->_capacity_values[initc + j]  = icap->second;
        j++;
      }
      j = 0;
      const std::map<double, double> *kappa_map = imat.second.retKappaMap();
      int initk = (*mt)->_kappa_inits[i];
      for(std::map<double,double>::const_iterator ikap = kappa_map->begin(); ikap != kappa_map->end(); ++ikap){
        (*mt)->_kappa_temps[initk + j] = ikap->first;
        (*mt)->_kappa_values[initk + j]  = ikap->second;
        j++;
      }
      i++;
    }
}
void freeMaterialTableMemory(MaterialTable **mt)
{
  free((*mt)->capacity);
  free((*mt)->kappa);
  free((*mt)->beta);
  free((*mt)->density);
  free((*mt)->_capacity_counters);
  free((*mt)->_capacity_inits);
  free((*mt)->_kappa_counters);
  free((*mt)->_kappa_inits);
  free((*mt)->_capacity_temps);
  free((*mt)->_capacity_values);
  free((*mt)->_kappa_temps);
  free((*mt)->_kappa_values);
}
