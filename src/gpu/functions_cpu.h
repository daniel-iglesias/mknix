#ifndef FUNCTIONS_CPU_H
#define FUNCTIONS_CPU_H

double computeCellAverageTemperature();

class ShapeFunctionTable{
public:
  int dim;
  double *phis; // i derivative odrer, j node
  int number_points;
  int* gp_ids;
};

double getPhi_from_table(ShapeFunctionTable &shapeFunctionTable,
                         int gp_id,
                         int derivative_order,
                         int node_id);

void setPhi_into_table(double phiValue,
                       ShapeFunctionTable &shapeFunctionTable,
                       int gp_id,
                       int derivative_order,
                       int node_id);

class CellTable{
public:
  int number_cells;
  int *material_id;

  double *alphas;
  int *gaussPoints_count;
  int *gaussPoint_init;
  int *gaussPoints_ids;
  /*int *gaussPoints_MC_count;
  int *gaussPoint_MC_init;
  int *gaussPoints_MC_ids;*/
  double *jacobians;
  /*int *bodyPoints_count;
  int *bodyPoint_init;
  int *bodyPoints_ids;*/
  double *dcs;
};

//class MaterialTable
class MaterialTable{
public:
  int number_materials;
  double *capacity;
  double *kappa;
  double *beta;
  double *density;
  double *m_capacity_temps;
  double *m_capacity_values;
  int *m_capacity_counters;
  int *m_capacity_inits;
  double *m_kappa_temps;
  double *m_kappa_values;
  int *m_kappa_counters;
  int *m_kappa_inits;
};

//Material functions
void addMaterialToTable();

double getMaterialKappa (MaterialTable &materials,
                         int material_id,
                         double average_temperature);

double getMaterialDensity (MaterialTable &materials,
                          int material_id);

double getMaterialCapacity (MaterialTable &materials,
                            int material_id,
                            double average_temperature);

#endif //FUNCTIONS_CPU_H
