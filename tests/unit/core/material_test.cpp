#include "core_test_support.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <core/material.h>
#include <core/point.h>

namespace
{

using core_test::expectNear;

class FixedTemperaturePoint : public mknix::Point
{
public:
    FixedTemperaturePoint(double temperature, double x, double y, double z)
        : Point(3, 0, x, y, z, 0.0, 0.0)
        , temperature_(temperature)
    {
    }

    double getTemp() const override
    {
        return temperature_;
    }

private:
    double temperature_;
};

void testThermalPropertiesAndInterpolation()
{
    mknix::Material material;
    material.setThermalProps(1000.0, 4.0, 0.001, 2.5);

    expectNear(material.getCapacity(), 1000.0, "constant heat capacity");
    expectNear(material.getKappa(), 4.0, "constant conductivity");
    expectNear(material.getBeta(), 0.001, "constant expansion coefficient");
    expectNear(material.getDensity(), 2.5, "constant density");

    material.addThermalCapacity(0.0, 800.0);
    material.addThermalCapacity(100.0, 1200.0);
    material.addThermalConductivity(0.0, 2.0);
    material.addThermalConductivity(100.0, 6.0);
    material.addThermalExpansion(0.0, 0.0005);
    material.addThermalExpansion(100.0, 0.0015);
    material.addThermalDensity(0.0, 3.0);
    material.addThermalDensity(100.0, 1.0);

    expectNear(material.getCapacity(50.0), 1000.0, "interpolated heat capacity");
    expectNear(material.getKappa(50.0), 4.0, "interpolated conductivity");
    expectNear(material.getBeta(50.0), 0.001, "interpolated expansion coefficient");
    expectNear(material.getDensity(50.0), 2.0, "interpolated density");
    expectNear(material.getCapacity(-10.0), 800.0, "lower-bound heat capacity");
    expectNear(material.getCapacity(110.0), 1200.0, "upper-bound heat capacity");
}

void testTwoDimensionalElasticMatrices()
{
    mknix::Material material;
    material.setMechanicalProps(2, 200.0, 0.25, 7.0);

    lmx::DenseMatrix<double>& linear = material.getD();
    lmx::DenseMatrix<double>& nonlinear = material.getC();
    const double lambda = 80.0;
    const double mu = 80.0;

    expectNear(linear.readElement(0, 0), 200.0 / (1.0 - 0.25 * 0.25), "2D D(0,0)");
    expectNear(linear.readElement(0, 1), 50.0 / (1.0 - 0.25 * 0.25), "2D D(0,1)");
    expectNear(linear.readElement(2, 2), 75.0 / (1.0 - 0.25 * 0.25), "2D D(2,2)");
    expectNear(nonlinear.readElement(0, 0), lambda + 2.0 * mu, "2D C(0,0)");
    expectNear(nonlinear.readElement(0, 1), lambda, "2D C(0,1)");

    int i = 0;
    int j = 1;
    int k = 0;
    int l = 1;
    expectNear(material.getCsym(i, j, k, l), mu, "minor-symmetric shear tensor component");
}

void testThreeDimensionalElasticMatrices()
{
    mknix::Material material;
    material.setMechanicalProps(3, 120.0, 0.2, 5.0);

    lmx::DenseMatrix<double>& linear = material.getD();
    lmx::DenseMatrix<double>& nonlinear = material.getC();
    const double factor = 120.0 / ((1.0 + 0.2) * (1.0 - 2.0 * 0.2));
    const double lambda = (0.2 * 120.0) / ((1.0 + 0.2) * (1.0 - 2.0 * 0.2));
    const double mu = 50.0;

    expectNear(linear.readElement(0, 0), factor * 0.8, "3D D(0,0)");
    expectNear(linear.readElement(0, 1), factor * 0.2, "3D D(0,1)");
    expectNear(linear.readElement(5, 5), factor * 0.3, "3D D(5,5)");
    expectNear(nonlinear.readElement(0, 0), lambda + 2.0 * mu, "3D C(0,0)");
    expectNear(nonlinear.readElement(2, 1), lambda, "3D C(2,1)");
    expectNear(nonlinear.readElement(4, 4), mu, "3D C(4,4)");
}

void testPorosityResistanceLoadAndHistory()
{
    mknix::Material material;
    material.setPorosityProps(9.0, 10.0, 4.0);
    material.addVariableResistance(0.0, 2.0);
    material.addVariableResistance(10.0, 4.0);

    FixedTemperaturePoint point(22.0, 5.0, 0.0, 0.0);
    expectNear(material.getPorosityResistance(&point), 3.0, "interpolated porosity resistance");
    expectNear(material.computePorosityLoad(&point), 4.0, "porosity heat load");
    expectNear(material.getPorosityFluidTemp(), 11.0, "capacity-adjusted fluid temperature");

    material.update(1);
    expectNear(material.getFluidTempHistory().back(), 11.0, "committed fluid temperature");
    expectNear(static_cast<double>(material.getFluidTempHistory().size()), 2.0,
               "fluid temperature history length");

    material.setPorosityFluidTemp(20.0);
    material.update(0);
    expectNear(material.getPorosityFluidTemp(), 10.0, "fluid temperature rollback");
}

} // namespace

void testMaterial()
{
    testThermalPropertiesAndInterpolation();
    testTwoDimensionalElasticMatrices();
    testThreeDimensionalElasticMatrices();
    testPorosityResistanceLoadAndHistory();
}