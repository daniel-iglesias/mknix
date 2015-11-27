//-- Licencia --
#include "LMX/lmx.h"
#include "cell.h"
#include "material.h"
#include "gausspoint.h"
#include "node.h"
#include "loadthermalbody.h"

namespace mknix {

Cell::Cell()
{
}


Cell::Cell( Material& material_in,
            std::string formulation_in,
            double alpha_in,
            int nGPoints_in
          )
    : mat(&material_in)
    , formulation(formulation_in)
    , alpha(alpha_in)
    , nGPoints(nGPoints_in)
    , dc(0)
{
}

Cell::~Cell()
{
  std::vector< GaussPoint* >::iterator it_gPoints;
  for(it_gPoints=gPoints.begin();
      it_gPoints!=gPoints.end();
      ++it_gPoints){
    delete(*it_gPoints);
  }
    for(it_gPoints=gPoints_MC.begin();
      it_gPoints!=gPoints_MC.end();
      ++it_gPoints){
    delete(*it_gPoints);
  }
}

// Only for Meshfree Cells, function is specialized for FEM elements
void Cell::initialize( std::vector<Node*> & nodes_in )
{
    // This function can be joined with assembleGaussPoints so the Gpoints are iterated only once...
    std::vector<GaussPoint*>::iterator it_gp;

    for ( it_gp = gPoints.begin();
            it_gp != gPoints.end();
            ++it_gp)
    {
        gPoints_MC.push_back(*it_gp); // use same GP for all matrices
        (*it_gp)->findSupportNodes( nodes_in );
    }
    // Set the dc and alpha parameteres for cell nodes. This way, the values will
    // be greater than zero only for meshfree nodes which need shapefunctions to be
    // calculated. The type of shape function is also set here.
//   std::vector<Node*>::iterator it_p;
//   for ( it_p = nodes.begin();
//         it_p != nodes.end();
//         ++it_p)
//   {
//     (*it_p)->setAlphai( alpha );
//     (*it_p)->setDc( dc );
//     if( formulation == "RPIM" )
//       (*it_p)->setShapeFunType( "RBF" );
//     else if( formulation == "EFG" )
//       (*it_p)->setShapeFunType( "MLS" );
//   }

}


void Cell::computeShapeFunctions(  )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        if( formulation == "RPIM" )
            (*it)->shapeFunSolve( "RBF", 1.03 );
        else if( formulation == "EFG" )
            (*it)->shapeFunSolve( "MLS", 1.03 );
    }
}


void Cell::computeCapacityGaussPoints(  )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints_MC.begin();
            it != gPoints_MC.end();
            ++it)
    {
        (*it)->computeCij( );
    }
}

void Cell::assembleCapacityGaussPoints( lmx::Matrix< data_type > & globalCapacity )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints_MC.begin();
            it != gPoints_MC.end();
            ++it)
    {
        (*it)->assembleCij( globalCapacity );
    }
}


void Cell::computeConductivityGaussPoints(  )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->computeHij( );
    }
}

void Cell::assembleConductivityGaussPoints( lmx::Matrix< data_type > & globalConductivity )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->assembleHij( globalConductivity );
    }
}


void Cell::computeQextGaussPoints( LoadThermalBody* loadThermalBody_in )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->computeQext( loadThermalBody_in );
    }
}

void Cell::assembleQextGaussPoints( lmx::Vector< data_type > & globalQext )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->assembleQext( globalQext );
    }
}


void Cell::computeMGaussPoints(  )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints_MC.begin();
            it != gPoints_MC.end();
            ++it)
    {
        (*it)->computeMij( );
    }
}


void Cell::assembleMGaussPoints( lmx::Matrix< data_type > & globalMass )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints_MC.begin();
            it != gPoints_MC.end();
            ++it)
    {
        (*it)->assembleMij( globalMass );
    }
}


void Cell::computeFintGaussPoints(  )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->computeFint( );
    }
}


void Cell::computeNLFintGaussPoints(  )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->computeNLFint( );
    }
}


void Cell::assembleFintGaussPoints( lmx::Vector< data_type > & globalFint )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->assembleFint( globalFint );
    }
}


void Cell::computeFextGaussPoints(  )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints_MC.begin();
            it != gPoints_MC.end();
            ++it)
    {
        (*it)->computeFext( );
    }
}


void Cell::assembleFextGaussPoints( lmx::Vector< data_type > & globalFext )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints_MC.begin();
            it != gPoints_MC.end();
            ++it)
    {
        (*it)->assembleFext( globalFext );
    }
}


void Cell::computeKGaussPoints(  )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->computeKij( );
    }
}


void Cell::computeNLKGaussPoints(  )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->computeNLKij( );
    }
}


void Cell::assembleKGaussPoints( lmx::Matrix< data_type > & globalTangent )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->assembleKij( globalTangent );
    }
}


void Cell::assembleRGaussPoints( lmx::Vector< data_type > & globalStress,
                                 int firstNode
                               )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->computeStress( );
        (*it)->assembleRi( globalStress, firstNode );
    }
}


void Cell::assembleNLRGaussPoints( lmx::Vector< data_type > & globalStress,
                                   int firstNode
                                 )
{
    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        (*it)->computeNLStress( );
        (*it)->assembleRi( globalStress, firstNode );
    }
}


double Cell::calcPotentialEGaussPoints( const lmx::Vector<data_type> & q )
{
    double potentialEnergy = 0;

    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        potentialEnergy += (*it)->calcPotentialE( q );
    }
    return potentialEnergy;
}


double Cell::calcKineticEGaussPoints( const lmx::Vector<data_type> & qdot )
{
    double kineticEnergy = 0;

    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        kineticEnergy += (*it)->calcKineticE( qdot );
    }
    return kineticEnergy;
}


double Cell::calcElasticEGaussPoints( )
{
    double elasticEnergy = 0;

    for ( std::vector<GaussPoint*>::iterator it = gPoints.begin();
            it != gPoints.end();
            ++it)
    {
        elasticEnergy += (*it)->calcElasticE( );
    }
    return elasticEnergy;
}


void Cell::outputConnectivityToFile(std::ofstream* outfile)
{
  std::vector< Point* >::iterator it_points;
  *outfile << "\t\t\t";
  for(it_points=bodyPoints.begin();
      it_points!=bodyPoints.end();
      ++it_points){
    *outfile << (*it_points)->getNumber() << " ";
  }
  *outfile << std::endl;
}


void Cell::gnuplotOutStress( std::ofstream & gptension )
{
    int counter;
    for(std::vector<GaussPoint*>::iterator it=gPoints.begin();
            it!=gPoints.end();
            ++it)
    {
        ++counter;
        (*it)->gnuplotOutStress( gptension );
        if (counter%4 == 0) gptension << endl;
    }
}

} //Namespace mknix
