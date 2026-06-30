//-- Licencia --
#include "gausspoint.h"
#include "node.h"
#include "shapefunctionRBF.h"
#include "shapefunctionMLS.h"

#include <simulation/simulation.h>
#include <system/system.h>
#include <system/loadthermalbody.h>

namespace mknix
{

/**
 * @brief Default constructor for GaussPoint.
 */
GaussPoint::GaussPoint()
{
}


/**
 * @brief Constructs a 2D Gauss point with coordinates (x, y).
 * @param dim_in Spatial dimension.
 * @param alpha_in Influence radius scaling factor.
 * @param weight_in Quadrature weight.
 * @param jacobian_in Jacobian of the mapping to the reference domain.
 * @param mat_in Pointer to the material at this point.
 * @param num_in Index of this Gauss point within its cell.
 * @param coor_x X coordinate.
 * @param coor_y Y coordinate.
 * @param dc_in Characteristic nodal spacing.
 * @param stressPoint_in True if this point is used for stress smoothing.
 */
GaussPoint::GaussPoint(int dim_in,
                       double alpha_in,
                       double weight_in,
                       double jacobian_in,
                       Material * mat_in,
                       int num_in,
                       double coor_x,
                       double coor_y,
                       double dc_in,
                       bool stressPoint_in
                      )
    : Point(dim_in, num_in, coor_x, coor_y, 0., jacobian_in, alpha_in, dc_in)
    , weight(weight_in)
    , mat(mat_in)
    , stressPoint(stressPoint_in)
{
}

/**
 * @brief Constructs a 3D Gauss point with coordinates (x, y, z).
 * @param dim_in Spatial dimension.
 * @param alpha_in Influence radius scaling factor.
 * @param weight_in Quadrature weight.
 * @param jacobian_in Jacobian of the mapping.
 * @param mat_in Pointer to the material.
 * @param num_in Index of this point within its cell.
 * @param coor_x X coordinate.
 * @param coor_y Y coordinate.
 * @param coor_z Z coordinate.
 * @param dc_in Characteristic nodal spacing.
 * @param stressPoint_in True if this point is used for stress smoothing.
 */
GaussPoint::GaussPoint(int dim_in,
                       double alpha_in,
                       double weight_in,
                       double jacobian_in,
                       Material * mat_in,
                       int num_in,
                       double coor_x,
                       double coor_y,
                       double coor_z,
                       double dc_in,
                       bool stressPoint_in
                      )
    : Point(dim_in, num_in, coor_x, coor_y, coor_z, jacobian_in, alpha_in, dc_in)
    , weight(weight_in)
    , mat(mat_in)
    , stressPoint(stressPoint_in)
{
}

/**
 * @brief Destructor for GaussPoint.
 */
GaussPoint::~GaussPoint()
{
}


/**
 * @brief Computes and stores shape function values at this point using the specified meshfree method.
 * @param type_in Shape function type: "RBF" for radial basis functions or "MLS" for moving least squares.
 * @param q_in Shape parameter (overridden internally to 0.5).
 */
void GaussPoint::shapeFunSolve(std::string type_in, double q_in)
{
    q_in = 0.5; // Original RBF
    if (!shapeFun)
    {
        if (type_in == "RBF")
        {
//             alphai = 3.5; // For the validation triangle, this works better: phi closer to 1.
            this->shapeFun = new ShapeFunctionRBF(supportNodesSize,
                                                  0,
                                                  0, // RBF type
                                                  alphai,
                                                  dc,
                                                  q_in,
                                                  this);
        }
        else if (type_in == "MLS")
        {
// 	    alphai=3.5; // For the validation triangle, this works better: phi closer to 1.
            this->shapeFun = new ShapeFunctionMLS(supportNodesSize,
                                                  1,
                                                  1, // weight type
                                                  alphai,
                                                  dc,
                                                  this);
        }
        // cout << "INFO AT shapeFunSolve IN GaussPoint: (x, y) = "
        //      << this->X << ", " << this->Y << endl;
        // cout << "\t alphai = " << alphai << ", "
        //      << "dc = " << dc << ", "
        //      << "q_in = " << q_in
        //      << endl;
        // cout << "\t Number of Support Nodes = " << supportNodesSize << endl;

        shapeFun->calc();
    }
}


/**
 * @brief Computes the local thermal capacity (heat capacity) matrix contribution C at this Gauss point.
 */
void GaussPoint::computeCij()
{
    // TODO: not sure why it's needed to be done here too, but for the moment is required for positive validation
    avgTemp = 0;
    for (auto i = 0u; i < supportNodesSize; ++i)
    {
        avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0, i);
    }
    double avgFactor = mat->getDensity(avgTemp) * mat->getCapacity(avgTemp) * weight * std::abs(jacobian);
    //////////////// Calculation of Capacity matrix:
    //////////////// M = rho * Cp * wg * N^T * N * |Jc|
    //////////////// Mij = rho * Cp * wg * Ni * Nj * |Jc| = M(i ,j)
    int j;
//   cout << mat->getDensity() << " " << mat->getCapacity() << " = density, capacity" << endl;
//     C.reset();
//    TODO: Select between lumped and consistent matrices as an input option
    for (auto i = 0u; i < supportNodesSize; ++i)
    {
        for (j=0; j<supportNodesSize; ++j) {
/////////////////////////////////
// Consistent Matrix option:   //
// Not reliable for coarse     //
// meshes, check for negative  //
// temperatures and refine     //
/////////////////////////////////
            C.writeElement( avgFactor * shapeFun->getPhi(0,i) * shapeFun->getPhi(0,j), i, j );
        }
/////////////////////////////////
// Do not use Lumped in ALICIA //
// It is also  buggy for 3D... //
/////////////////////////////////
// Lumped matrix (bug?):
        //  C.addElement( mat->getDensity() * mat->getCapacity() * weight * shapeFun->getPhi(0,i)
        //                     * shapeFun->getPhi(0,j) * std::abs(jacobian), i, i );
        // }
// Faster lumped matrix:
        // C.writeElement(
        //     mat->getDensity() * mat->getCapacity(supportNodes[i]->getTemp()) * weight * shapeFun->getPhi(0, i)
        //     * std::abs(jacobian), i, i);
// 	  cout << i << "," << j << " = "
// 	       << mat->getDensity() << "*"
// 	       << mat->getCapacity() << "*"
// 	       << weight  << "*"
// 	       << shapeFun->getPhi(0,i) << "*"
// 	       << shapeFun->getPhi(0,j)  << "*"
// 	       << std::abs(jacobian)  << " = "
// 	       << C.readElement(i,j) << endl;
    }
//     HRZ lumped matrix [Hinton et al. (1976)] (if implemented, should be in Cell class)

}

/**
 * @brief Computes the local thermal conductivity matrix contribution H at this Gauss point.
 */
void GaussPoint::computeHij()
{
    int max_deriv_index = dim + 1;
    H.reset();
    avgTemp = 0;
    for (auto i = 0u; i < supportNodesSize; ++i)
    {
        avgTemp += supportNodes[i]->getTemp() * shapeFun->getPhi(0, i);
    }
    double avgFactor = mat->getKappa(avgTemp) * weight * std::abs(jacobian);
    // Hij = wg * grad(N_j) * kappa * grad(N_I) * |Jc|
    for (auto i = 0u; i < supportNodesSize; ++i)
    {
        for (auto j = 0u; j < supportNodesSize; ++j)
        {
// 	  if(supportNodes[i]->getThermalNumber() == 39){
// 	    cout << "KAPPA in " << i << "," << j << " = " ;
// 	    cout << 0.5*( mat->getKappa(supportNodes[i]->getTemp()) + mat->getKappa(supportNodes[j]->getTemp()) );
// 	    cout << " with nodes temperatures:  ";
// 	    cout << supportNodes[i]->getTemp() << " and ";
// 	    cout << supportNodes[j]->getTemp() << endl;
// 	  }
            for (auto m = 1; m < max_deriv_index; ++m)
            {
//         for ( n=1; n<max_deriv_index; ++n ){
                //           K(2*i + m, 2*j + n) = Kij(m,n);
                H.addElement( (shapeFun->getPhi(m, i)
                               // 					 * mat->getKappa(supportNodes[i]->getTemp())
                               // 					 + .5*mat->getKappa(supportNodes[j]->getTemp())/*.readElement(m,n)*/
//                                        * mat->getKappa(avgTemp) /*.readElement(m,n)*/
                               // 					 * 0.5*( mat->getKappa(supportNodes[i]->getTemp()) + mat->getKappa(supportNodes[j]->getTemp()) )/*.readElement(m,n)*/
                               * shapeFun->getPhi(m, j)) * avgFactor,
                              i,
                              j);
// 	  cout << i << "," << j << " = "
// 	       << mat->getDensity() << "*"
// 	       << mat->getKappa() << "*"
// 	       << weight  << "*"
// 	       << shapeFun->getPhi(m,i) << "*"
// 	       << shapeFun->getPhi(m,j)  << "*"
// 	       << jacobian  << " = "
// 	       << H.readElement(i,j) << endl;
//         }
            }
//         H.writeElement( weight * mat->getKappa() * jacobian *
//                         (  shapeFun->getPhi(1,i) * shapeFun->getPhi(1,j) +
//                          + shapeFun->getPhi(2,i) * shapeFun->getPhi(2,j)
// 			),
//                                   i,
//                                   j);
        }
    }
//     cout << "GP (" << this->getX() << ", " << this->getY() << ")" << endl;
//     cout << "H = " << H << endl;
}

// Todo: add somewhere the coordinate for interpolation of heat load and porosity resistance. It now uses only x.
// Also, it would be ideal to divide this function into two parts so that we separate 
/**
 * @brief Computes the external thermal load vector contribution Qext at this Gauss point,
 *        including volumetric heat sources and porosity cooling.
 * @param loadThermalBody_in Vector of body thermal load objects.
 */
void GaussPoint::computeQext(std::vector<LoadThermalBody*> loadThermalBody_in)
{
    double load;
    // Only execute if the load vector is not empty or if the porosity resistance is greater than zero. 
    if ( !loadThermalBody_in.empty() || mat->isPorous() ) {
        // Compute the coordinate for the external volumetric power
        for (auto i = 0u; i < supportNodesSize; ++i)
        {
            // Reset the load
            load = 0.0;
            // Qi = wg * ( r- (T - Tl)/R ) * N_I * |Jc| , 
            // r-> external volumetric power, 
            // Tl->fluid temp, 
            // R->volumetric ratio
            for (auto& v_load : loadThermalBody_in){
                load += v_load->getLoadThermalBody(supportNodes[i]) ;
            }
            if(mat->isPorous()){
               load -= mat->computePorosityLoad(supportNodes[i]);
            }
            Qext.writeElement( weight * shapeFun->getPhi(0, i) * load * std::abs(jacobian), i );
        }
    }
}


/**
 * @brief Assembles the local capacity matrix C into the global thermal capacity matrix.
 * @param globalCapacity Global thermal capacity matrix to be updated.
 */
void GaussPoint::assembleCij(lmx::Matrix<data_type>& globalCapacity)
{
    for (auto i = 0u; i < supportNodesSize; ++i)
    {
        for (auto j = 0u; j < supportNodesSize; ++j)
        {
            globalCapacity.addElement(C.readElement(i, j),
                                      supportNodes[i]->getThermalNumber(),
                                      supportNodes[j]->getThermalNumber()
                                     );
        }
    }
//    cout << globalCapacity << endl;
}

/**
 * @brief Assembles the local conductivity matrix H into the global conductivity matrix.
 * @param globalConductivity Global conductivity matrix to be updated.
 */
void GaussPoint::assembleHij(lmx::Matrix<data_type>& globalConductivity)
{
    for (auto i = 0u; i < supportNodesSize; ++i)
    {
        for (auto j = 0u; j < supportNodesSize; ++j)
        {
            globalConductivity.addElement(H.readElement(i, j),
                                          supportNodes[i]->getThermalNumber(),
                                          supportNodes[j]->getThermalNumber()
                                         );
        }
    }
//    cout << globalConductivity << endl;
}

/**
 * @brief Assembles the local external heat vector Qext into the global heat load vector.
 * @param globalHeat Global external heat load vector to be updated.
 */
void GaussPoint::assembleQext(lmx::Vector<data_type>& globalHeat)
{
    for (auto i = 0u; i < supportNodesSize; ++i)
    {
        globalHeat.addElement(Qext.readElement(i),
                              supportNodes[i]->getThermalNumber()
                             );
    }
}


/**
 * @brief Outputs the Gauss point position and stress value to a gnuplot-compatible file.
 * @param gptension Output file stream for stress data.
 */
void GaussPoint::gnuplotOutStress(std::ofstream& gptension)
{
    gptension << X << " " << Y << " " << tension(0) << endl;
}


} //Namespace mknix
