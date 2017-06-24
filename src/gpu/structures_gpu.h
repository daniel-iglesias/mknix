
#ifndef STRUCTURES_GPU_H
#define STRUCTURES_GPU_H
#include "cpu/structures.h"

bool allocateTransferShapeFunctionTable(ShapeFunctionTable** device_shape_table,
                                        ShapeFunctionTable* host_shape_table);

bool freeDeviceShapeFunctionTable(ShapeFunctionTable** device_shape_table);

bool allocateTransferMaterialTable(MaterialTable** device_material_table,
                                   MaterialTable* host_material_table);

bool freeDeviceMaterialTable(MaterialTable** device_material_table);

#endif //STRUCTURES_GPU_H
