#ifndef FUNCTIONS_CPU_H
#define FUNCTIONS_CPU_H

#include "core/material.h"

double computeCellAverageTemperature();

struct ShapeFunctionTable{
  double *phis; // i derivative order, j node
  int support_node_size;
  int number_points;
};

double getPhi_from_table(ShapeFunctionTable *shapeFunctionTable,
                         int point_id,
                         int local_node_id);

void setPhi_into_table(double phiValue,
                       ShapeFunctionTable *shapeFunctionTable,
                       int point_id,
                       int local_node_id);

//shapefunctions functions_cpuvoid
void init_shape_functions_table(ShapeFunctionTable **shapeFunctionTable,
                                int support_node_size,
                                int number_points);

/*struct CellTable{
  int number_cells;
  int *material_id;

  double *alphas;
  int *gaussPoints_count;
  int *gaussPoint_init;
  int *gaussPoints_ids;
  double *jacobians;
  double *dcs;
};*/


//class MaterialTable
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


//Material functions
void addMaterialToTable();

double getMaterialKappa (MaterialTable *materials,
                         int material_id,
                         double average_temperature);

double getMaterialDensity (MaterialTable *materials,
                          int material_id);

double getMaterialCapacity (MaterialTable *materials,
                            int material_id,
                            double average_temperature);
//
void setupMaterialTables(MaterialTable **mt,
                        std::map<int, mknix::Material> &materials);

#endif //FUNCTIONS_CPU_H
