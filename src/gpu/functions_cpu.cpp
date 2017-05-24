//some helper functions, arrangeds to fulfill Structure Of Arrays (SOA) philoshophy
#include "functions_cpu.h"
#include <stdlib.h>
#include <map>

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
    double this_val = reference_values[init_position + i];
    if(query_value <= this_val){//bounded here
      upper_index = i;
      break;
    }
  }
  if(upper_index == 0) return sought_values[init_position];
  upper_bounded = true;
  if(!upper_bounded) return sought_values[init_position + counter];//not bound found so return last
  //so we have a bound
  int lower_index = upper_index - 1;
  double delta = (query_value - reference_values[init_position + lower_index]) / (reference_values[init_position + upper_index] - reference_values[init_position + lower_index]);
  return delta * sought_values[init_position + upper_index] + (1.0 - delta) * sought_values[init_position + lower_index];
}

/*
struct MaterialTable{
  int number_materials;
  double *capacity;
  double *kappa;
  double *beta;
  double *density;
  double *_capacity_temps;
  double *_capacity_values;
  int *_capacity_counters;
  int *_capacity_inits;
  double *_kappa_temps;
  double *_kappa_values;
  int *_kappa_counters;
  int *_kappa_inits;
};

*/

void debug_printMaterialTable(MaterialTable *materials)
{
  std::cout<< std::endl << "--------------- DEBUG PRINTING MATERIAL TABLES ---------------" << std::endl;
  int nmat = materials->number_materials;
  std::cout << "Total Number of Materials = " << nmat << std::endl;
  for(int imat = 0; imat < nmat; imat++)
  {
    std::cout << "--------- Material " << imat <<"---------" <<  std::endl;
    std::cout << "capacity = " << materials->capacity[imat] << std::endl;
    std::cout << "kappa = " << materials->kappa[imat] << std::endl;
    std::cout << "beta = " << materials->beta[imat] << std::endl;
    std::cout << "density = " << materials->density[imat] << std::endl;
    int cap_pairs = materials->_capacity_counters[imat];
    int cap_init = materials->_capacity_inits[imat];
    std::cout << "capacity temps and vals pairs: "<< std::endl;
    for (int icap = 0; icap < cap_pairs; icap++){
      std::cout << materials->_capacity_temps[cap_init + icap] << "/" << materials->_capacity_values[cap_init + icap]<< "; ";
    }
    double check_cap = interpolate1D(375,
                                     materials->_capacity_temps,
                                     materials->_capacity_values,
                                     cap_init,
                                     cap_pairs);
    std::cout << "Checking interpolation with " << 375 << " results in " << check_cap << std::endl;
    int kap_pairs = materials->_kappa_counters[imat];
    int kap_init = materials->_kappa_inits[imat];
    std::cout << std::endl << "kappa temps and vals pairs: "<< std::endl;
    for (int ikap = 0; ikap < kap_pairs; ikap++){
      std::cout << materials->_kappa_temps[kap_init + ikap] << "/" << materials->_kappa_values[kap_init + ikap]<< "; ";
    }
    double check_kap = interpolate1D(355,
                                     materials->_kappa_temps,
                                     materials->_kappa_values,
                                     kap_init,
                                     kap_pairs);
    std::cout << "Checking interpolation with " << 355 << " results in " << check_kap << std::endl;

    std::cout << std::endl << "--------- End of Material " << imat <<"---------" <<  std::endl;
  }


  std::cout<< std::endl << "------------ END OF DEBUG PRINTING MATERIAL TABLES ------------" << std::endl;

}

double getPhi_from_table(ShapeFunctionTable *shapeFunctionTable,
                         int point_id,
                         int dim,
                         int local_node_id)
{
  int support_node_size = shapeFunctionTable->support_node_size;
  int max_deriv =  shapeFunctionTable->number_derivatives;
  return shapeFunctionTable->phis[point_id * support_node_size * max_deriv + dim * support_node_size + local_node_id];
}

void setPhi_into_table(double phiValue,
                       ShapeFunctionTable *shapeFunctionTable,
                       int point_id,
                       int dim,
                       int local_node_id)
{
  int support_node_size = shapeFunctionTable->support_node_size;
  int max_deriv =  shapeFunctionTable->number_derivatives;
  shapeFunctionTable->phis[point_id * support_node_size * max_deriv + dim * support_node_size + local_node_id] = phiValue;
}

void init_shape_functions_table(ShapeFunctionTable **shapeFunctionTable,
                                int support_node_size,
                                int number_derivatives,
                                int number_points)
{
  (*shapeFunctionTable)->support_node_size = support_node_size;
  (*shapeFunctionTable)->number_points = number_points;
  (*shapeFunctionTable)->phis = (double*) malloc(number_points * number_derivatives * support_node_size *  sizeof(double));
}

void free_shape_functions_table(ShapeFunctionTable **shapeFunctionTable)
{
  free((*shapeFunctionTable)->phis);
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
    /*i = 0;
    for(auto const &imat : materials){
      (*mt)->capacity[i] = imat.second.retCapacity();
      (*mt)->_capacity_counters[i]= imat.second.retCapacityMapSize();
      (*mt)->kappa[i] = imat.second.retKappa();
      (*mt)->_kappa_counters[i]= imat.second.retKappaMapSize();
      i++;
    }*/
    (*mt)->_capacity_inits[0] = 0;
    (*mt)->_kappa_inits[0] = 0;
    int capacity_tot_vals = (*mt)->_capacity_counters[0];
    int kappa_tot_vals = (*mt)->_kappa_counters[0];
    for(int imat = 1; imat < nMat; imat++){//prefix sum scan
      capacity_tot_vals += (*mt)->_capacity_counters[imat];
      kappa_tot_vals += (*mt)->_kappa_counters[imat];
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

//class MaterialTable
/*struct ThermalBoundaryTable{
  int number_thermal_boundaries;
  double *_load_coord;
  double *_load_values;
  int *_load_counters;
  int *_load_inits;
  double *_time_elapsed;
  double *_time_values;
  int *_time_counters;
  int *_time_inits;
};*/

void debug_printThermalBoundaryTable(ThermalBoundaryTable *thermalBoundaries)
{

}

double getThermalBoundaryLoad (ThermalBoundaryTable *thermalBoundaries,
                               int thermal_boundary_id,
                               double coord_x)
{
  int n_vals =  thermalBoundaries->_load_counters[thermal_boundary_id];
  int init_vals =  thermalBoundaries->_load_inits[thermal_boundary_id];
  return interpolate1D(coord_x,
                       thermalBoundaries->_load_coord,
                       thermalBoundaries->_load_values,
                       init_vals,
                       n_vals);
}

double getThermalBoundaryTime (ThermalBoundaryTable *thermalBoundaries,
                               int thermal_boundary_id,
                               double simulation_time)
{
  int n_vals =  thermalBoundaries->_time_counters[thermal_boundary_id];
  int init_vals =  thermalBoundaries->_time_inits[thermal_boundary_id];
  return interpolate1D(simulation_time,
                       thermalBoundaries->_time_elapsed,
                       thermalBoundaries->_time_values,
                       init_vals,
                       n_vals);
}

void setupThermalBoundaryTable(ThermalBoundaryTable **tbs,
                               std::map<int, mknix::LoadThermalBoundary1D> &boundaries)
{
  std::cout << "setting up SOA thermalBoundaries tables" << std::endl;
  int nBound = boundaries.size();
  std::cout << "Allocating memory for "<< nBound << " boundaries" << std::endl;
  (*tbs)->number_thermal_boundaries = nBound;
  (*tbs)->_load_counters = (int*) malloc(nBound * sizeof(int));
  (*tbs)->_load_inits = (int*) malloc(nBound * sizeof(int));
  (*tbs)->_time_counters = (int*) malloc(nBound * sizeof(int));
  (*tbs)->_time_inits = (int*) malloc(nBound * sizeof(int));
  int i = 0;
  std::cout << "Copy variables for materials" << std::endl;
  for(auto  &ibound : boundaries){
    (*tbs)->_load_counters[i]= ibound.second.getLoadMapSize();
    (*tbs)->_time_counters[i]= ibound.second.getTimeMapSize();
    i++;
  }
  (*tbs)->_load_inits[0] = 0;
  (*tbs)->_time_inits[0] = 0;
  int load_tot_vals = (*tbs)->_load_counters[0];
  int time_tot_vals = (*tbs)->_time_counters[0];
  for(int ib = 1; ib < nBound; ib++){//prefix sum scan
    load_tot_vals += (*tbs)->_load_counters[ib];
    time_tot_vals += (*tbs)->_time_counters[ib];
    (*tbs)->_load_inits[ib] = (*tbs)->_load_inits[ib - 1] + (*tbs)->_load_counters[ib - 1];
    (*tbs)->_time_inits[ib] = (*tbs)->_time_inits[ib - 1] + (*tbs)->_time_counters[ib - 1];
  }
  (*tbs)->_load_coord = (double*) malloc(load_tot_vals * sizeof(double));
  (*tbs)->_load_values = (double*) malloc(load_tot_vals * sizeof(double));
  (*tbs)->_time_elapsed = (double*) malloc(time_tot_vals * sizeof(double));
  (*tbs)->_time_values = (double*) malloc(time_tot_vals* sizeof(double));
  //now we copy all values
  i = 0;
  for(auto  &ibound : boundaries){
    const std::map<double, double> *load_map = ibound.second.getLoadMap();
    int initl = (*tbs)->_load_inits[i];
    //for(int icap = 0;  icap < cap_map->size(); icap++){
    int j =0;
    //for (auto const& icap:cap_map){
    for(std::map<double,double>::const_iterator ilo = load_map->begin(); ilo != load_map->end(); ++ilo){
      (*tbs)->_load_coord[initl + j] = ilo->first;
      (*tbs)->_load_values[initl + j]  = ilo->second;
      j++;
    }
    j = 0;
    const std::map<double, double> *time_map = ibound.second.getTimeMap();
    int initt = (*tbs)->_time_inits[i];
    for(std::map<double,double>::const_iterator iti = time_map->begin(); iti != time_map->end(); ++iti){
      (*tbs)->_time_elapsed[initt + j] = iti->first;
      (*tbs)->_time_values[initt + j]  = iti->second;
      j++;
    }
    i++;
  }

}
void freeThermalBoundaryTableMemory(ThermalBoundaryTable **tbs)
{
  free((*tbs)->_load_coord);
  free((*tbs)->_load_values);
  free((*tbs)->_load_counters);
  free((*tbs)->_load_inits);
  free((*tbs)->_time_elapsed);
  free((*tbs)->_time_values);
  free((*tbs)->_time_counters);
  free((*tbs)->_time_inits);
}
