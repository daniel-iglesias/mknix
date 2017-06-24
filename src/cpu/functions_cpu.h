
#ifndef FUNCTIONS_CPU_H
#define FUNCTIONS_CPU_H

#include "core/material.h"
#include "system/loadthermalboundary1D.h"
#include "structures.h"

/**
 * SoA Structure containing shapefunction phi values
 * @param  {[type]} double* array             array of local phis
 * @param  {[type]} int size                  Number of nodes per gausspoints
 * @param  {[type]} int size                  Number of gausspoints
 */
/*struct ShapeFunctionTable{
  double *phis; // i derivative order, j node
  int support_node_size;
  int number_points;
  int number_derivatives; // dim + 1
};*/
/**
 * returns the shape function phi value for a given gausspoint and node
 * @param  {[type]} ShapeFunctionTable* struct  structure containing shapefunction phis
 * @param  {[type]} int id                       Gausspoint id
 * @param  {[type]} int id                       local node id
 */
double getPhi_from_table(ShapeFunctionTable *shapeFunctionTable,
                         int point_id,
                         int dim,
                         int local_node_id);
//
/**
 * Sets the shape function phi value for a given gausspoint and node
 * @param  {[type]} double value                value to set as phi
 * @param  {[type]} ShapeFunctionTable* struct  structure containing shapefunction phis
 * @param  {[type]} int id                       Gausspoint id
 * @param  {[type]} int id                       local node id
 */
void setPhi_into_table(double phiValue,
                       ShapeFunctionTable *shapeFunctionTable,
                       int point_id,
                       int deriv,
                       int local_node_id);

/**
 * Allocates memory for the shape function SoA struct
 * @param  {[type]} ShapeFunctionTable* struct  structure containing shapefunction phis
 * @param  {[type]} int size                     Number of support nodes per gausspoint
 * @param  {[type]} int id                       Number of gausspoints
 */
void init_shape_functions_table(ShapeFunctionTable **shapeFunctionTable,
                                int support_node_size,
                                int number_derivatives,
                                int number_points);
//
/**
 * Frees memory for the shape function SoA struct
 * @param  {[type]} ShapeFunctionTable* struct  structure containing shapefunction phis
 */
void free_shape_functions_table(ShapeFunctionTable **shapeFunctionTable);

/**
 * SoA Structure containing all materials data for coalescent access
 * @param  {[type]} int size                  Number of different materials in the structc
 * @param  {[type]} double* array             Array of capacity values for the materials
 * @param  {[type]} double* array             Array of kappa values for the materials
 * @param  {[type]} double* array             Array of beta values for the materials
 * @param  {[type]} double* array             Array of density values for the materials
 * @param  {[type]} double* array             Sparse matrix of temperatures map for capacity
 * @param  {[type]} double* array             Sparse matrix of values for capacity
 * @param  {[type]} int* array                Array of indices in sparse matrix capacity
 * @param  {[type]} int* array                Array of row counters in sparse matrix capacity
 * @param  {[type]} double* array             Sparse matrix of temperatures map for kappa
 * @param  {[type]} double* array             Sparse matrix of values for kappa
 * @param  {[type]} int* array                Array of indices in sparse matrix kappa
 * @param  {[type]} int* array                Array of row counters in sparse matrix kappa
 * @param  {[type]} int* array                Sparse matrix of values for capacity
 */
//class MaterialTable
/*struct MaterialTable{
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
};*/

void debug_printMaterialTable(MaterialTable *materials);

//Material functions
/**
 * returns kappa for a material at a given temperature
 * @param  {[type]} MaterialTable* struct  Structure containing materials data
 * @param  {[type]} int id                 Material id
 * @param  {[type]} double value           Temperature
 */
double getMaterialKappa (MaterialTable *materials,
                         int material_id,
                         double average_temperature);
//
/**
* returns density for a material
* @param  {[type]} MaterialTable* struct  Structure containing materials data
* @param  {[type]} int id                 Material id
 */
double getMaterialDensity (MaterialTable *materials,
                          int material_id);
//
/**
* returns capacity for a material at a given temperature
* @param  {[type]} MaterialTable* struct  Structure containing materials data
* @param  {[type]} int id                 Material id
* @param  {[type]} double value           Temperature
 */
double getMaterialCapacity (MaterialTable *materials,
                            int material_id,
                            double average_temperature);
//
/**
 * Allocates memory for the Materials SoA struct
 * @param  {[type]} MaterialTable* struct  structure containing materials data
 * @param  {[type]} std::map<int, mknix::Material> map Materials Objects map
 */
void setupMaterialTables(MaterialTable **mt,
                        std::map<int, mknix::Material> &materials);
//
/**
 * Frees memory for the materials SoA struct
 * @param  {[type]} MaterialTable* struct  structure containing materials data
 */
void freeMaterialTableMemory(MaterialTable **mt);
//

/**
 * SoA Structure containing all 1D Thermal Boundaries for coalescent access
 * @param  {[type]} int size                  Number of different thermal boundaries in the structc
 * @param  {[type]} double* array             Sparse matrix of coord map for load
 * @param  {[type]} double* array             Sparse matrix of values for load
 * @param  {[type]} int* array                Array of indices in sparse matrix load
 * @param  {[type]} int* array                Array of row counters in sparse matrix load
 * @param  {[type]} double* array             Sparse matrix of time map for time
 * @param  {[type]} double* array             Sparse matrix of values for time
 * @param  {[type]} int* array                Array of indices in sparse matrix time
 * @param  {[type]} int* array                Array of row counters in sparse matrix time
 */
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

void debug_printThermalBoundaryTable(ThermalBoundaryTable *thermalBoundaries);

//Material functions
/**
 * returns load for a 1D boundary at a given coord
 * @param  {[type]} ThermalBoundaryTable* struct  Structure containing materials data
 * @param  {[type]} int id                 Material id
 * @param  {[type]} double value           Temperature
 */
 double getThermalBoundaryLoad (ThermalBoundaryTable *thermalBoundaries,
                                int thermal_boundary_id,
                                double coord_x);

//
/**
* returns time factor for a 1D ThermalBoundary at a given time
* @param  {[type]} ThermalBoundaryTable* struct  Structure containing materials data
* @param  {[type]} int id                 Material id
* @param  {[type]} double value           Temperature
 */
 double getThermalBoundaryTime (ThermalBoundaryTable *thermalBoundaries,
                                int material_id,
                                double simulation_time);
//
/**
 * Allocates memory for the ThermalBoundary SoA struct
 * @param  {[type]} ThermalBoundaryTable** struct  structure containing boundaries data
 * @param  {[type]} std::map<int, mknix::LoadThermalBoundary1D> map 1D boundary Objects map
 */
void setupThermalBoundaryTable(ThermalBoundaryTable **thermalBoundaries,
                              std::map<int, mknix::LoadThermalBoundary1D*> boundaries);
//
/**
 * Frees memory for the ThermalBoundary SoA struct
 * @param  {[type]} ThermalBoundaryTable** struct  structure containing materials data
 */
void freeThermalBoundaryTableMemory(ThermalBoundaryTable **thermalBoundaries);
//



#endif //FUNCTIONS_CPU_H
