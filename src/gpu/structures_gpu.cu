#include "structures_gpu.h"
#include "cuda_helper.h"
/*struct ShapeFunctionTable{
  double *phis; // i derivative order, j node
  int support_node_size;
  int number_points;
  int number_derivatives; // dim + 1
};*/
bool allocateTransferShapeFunctionTable(ShapeFunctionTable** device_shape_table,
                                        ShapeFunctionTable* host_shape_table)
{
  int number_points      = host_shape_table->number_points;
  int number_derivatives = host_shape_table->number_derivatives;
  int support_node_size  = host_shape_table->support_node_size;

  CudaSafeCall(cudaMalloc((void**) device_shape_table, sizeof(ShapeFunctionTable)));

  CudaSafeCall(cudaMalloc((void**) &((*device_shape_table)->phis), number_points * number_derivatives * support_node_size * sizeof(double)));
//transfer
  CudaSafeCall(cudaMemcpy(&((*device_shape_table)->support_node_size),  &(host_shape_table->support_node_size),  1 * sizeof(int), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy(&((*device_shape_table)->number_points),      &(host_shape_table->number_points),      1 * sizeof(int), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy(&((*device_shape_table)->number_derivatives), &(host_shape_table->number_derivatives), 1 * sizeof(int), cudaMemcpyHostToDevice));

  CudaSafeCall(cudaMemcpy((*device_shape_table)->phis, host_shape_table->phis, number_points * number_derivatives * support_node_size * sizeof(double), cudaMemcpyHostToDevice));

  return true;
}

bool freeDeviceShapeFunctionTable(ShapeFunctionTable** device_shape_table)
{
  CudaSafeCall(cudaFree(&((*device_shape_table)->phis)));
  CudaSafeCall(cudaFree(*device_shape_table));
  return true;
}
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

bool allocateTransferMaterialTable(MaterialTable** device_material_table,
                                   MaterialTable* host_material_table)
{
  //first allocate teh Material Table in the GPU
  std::cout << "Allocate device Material Table" << std::endl;
  //CudaSafeCall(cudaMalloc((void**)&(*device_material_table->number_materials), sizeof(int)));
  int nMat = host_material_table->number_materials;
  std::cout << "Number materials " << nMat << std::endl;
  std::cout << "Allocate device Material Table param arrays" << std::endl;
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->capacity), nMat * sizeof(double)));
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->kappa),    nMat * sizeof(double)));
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->beta),     nMat * sizeof(double)));
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->density),  nMat * sizeof(double)));
 std::cout << "Allocate device Material Table counter and init arrays" << std::endl;
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->_capacity_counters), nMat * sizeof(int)));
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->_capacity_inits),    nMat * sizeof(int)));
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->_kappa_counters),    nMat * sizeof(int)));
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->_kappa_inits),     nMat * sizeof(int)));

  int capacity_tot_vals = 0;
  int kappa_tot_vals = 0;
  for(int imat = 0; imat < nMat; imat++){//prefix sum scan
      capacity_tot_vals += host_material_table->_capacity_counters[imat];
      kappa_tot_vals += host_material_table->_kappa_counters[imat];
  }
  std::cout << "kappa_tot_vals " << kappa_tot_vals << std::endl;
  std::cout << "Allocate device Material Table values and temps arrays" << std::endl;
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->_capacity_temps),  capacity_tot_vals * sizeof(double)));
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->_capacity_values), capacity_tot_vals * sizeof(double)));
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->_kappa_temps),     kappa_tot_vals * sizeof(double)));
  CudaSafeCall(cudaMalloc((void**) &((*device_material_table)->_kappa_values),    kappa_tot_vals * sizeof(double)));

  //now we transfer all memory
  std::cout << "Transfering all Material Table memory" << std::endl;
  //CudaSafeCall(cudaMemcpy(&((*device_material_table)->number_materials), &(host_material_table->number_materials),  1 * sizeof(int), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->capacity, host_material_table->capacity,  nMat * sizeof(double), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->kappa,    host_material_table->kappa,     nMat * sizeof(double), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->beta,     host_material_table->beta,      nMat * sizeof(double), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->density,  host_material_table->density,   nMat * sizeof(double), cudaMemcpyHostToDevice));

  CudaSafeCall(cudaMemcpy((*device_material_table)->_capacity_counters, host_material_table->_capacity_counters,  nMat * sizeof(int), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->_capacity_inits,    host_material_table->_capacity_inits,     nMat * sizeof(int), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->_kappa_counters,    host_material_table->_kappa_counters,     nMat * sizeof(int), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->_kappa_inits,       host_material_table->_kappa_inits,        nMat * sizeof(int), cudaMemcpyHostToDevice));

  CudaSafeCall(cudaMemcpy((*device_material_table)->_capacity_temps,  host_material_table->_capacity_temps,  capacity_tot_vals * sizeof(double), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->_capacity_values, host_material_table->_capacity_values, capacity_tot_vals * sizeof(double), cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->_kappa_temps,     host_material_table->_kappa_temps,     kappa_tot_vals * sizeof(double),    cudaMemcpyHostToDevice));
  CudaSafeCall(cudaMemcpy((*device_material_table)->_kappa_values,    host_material_table->_kappa_values,    kappa_tot_vals * sizeof(double),    cudaMemcpyHostToDevice));

  return true;
}

bool freeDeviceMaterialTable(MaterialTable** device_material_table)
{
  //first allocate teh Material Table in the GPU
  CudaSafeCall(cudaFree((*device_material_table)->capacity));
  CudaSafeCall(cudaFree((*device_material_table)->kappa));
  CudaSafeCall(cudaFree((*device_material_table)->beta));
  CudaSafeCall(cudaFree((*device_material_table)->density));

  CudaSafeCall(cudaFree((*device_material_table)->_capacity_counters));
  CudaSafeCall(cudaFree((*device_material_table)->_capacity_inits));
  CudaSafeCall(cudaFree((*device_material_table)->_kappa_counters));
  CudaSafeCall(cudaFree((*device_material_table)->_kappa_inits));

  CudaSafeCall(cudaFree((*device_material_table)->_capacity_temps));
  CudaSafeCall(cudaFree((*device_material_table)->_capacity_values));
  CudaSafeCall(cudaFree((*device_material_table)->_kappa_temps));
  CudaSafeCall(cudaFree((*device_material_table)->_kappa_values));

  CudaSafeCall(cudaFree(*device_material_table));

  return true;
}
/*
//class MaterialTable
struct ThermalBoundaryTable{
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
