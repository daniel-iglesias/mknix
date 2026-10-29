#include "core_test_support.h"

#include <core/node.h>
#include <core/point.h>
#include <core/shapefunction.h>
#include <core/shapefunctionlinear.h>
#include <core/shapefunctionlinear-x.h>
#include <core/shapefunctionMLS.h>
#include <core/shapefunctionMLS2D.h>
#include <core/shapefunctionMLS3D.h>
#include <core/shapefunctionRBF.h>
#include <core/shapefunctiontetrahedron.h>
#include <core/shapefunctiontriangle.h>
#include <core/shapefunctiontriangle3D.h>

namespace
{

class TestShapeFunction : public mknix::ShapeFunction
{
public:
    TestShapeFunction()
    {
        phi.resize(1, 1);
    }

    void calc() override
    {
    }
};

void testShapeFunctionBase()
{
    TestShapeFunction shape;
    shape.setPhi(0.75, 0, 0);
    core_test::expectNear(shape.getPhi(0, 0), 0.75, "ShapeFunction stores and returns a value");
}

void addNodes(mknix::Point& point, std::vector<mknix::Node*>& nodes)
{
    for (std::vector<mknix::Node*>::iterator node = nodes.begin(); node != nodes.end(); ++node)
    {
        point.addSupportNode(*node);
    }
}

void testLinearShapeFunctions()
{
    mknix::Node start(0, 0.0, 0.0, 0.0);
    mknix::Node end(1, 2.0, 0.0, 0.0);
    std::vector<mknix::Node*> nodes = {&start, &end};

    mknix::Point point1D(1, 2, 0.5, 0.0, 0.0, 1.0, 2.0);
    addNodes(point1D, nodes);
    mknix::ShapeFunctionLinear linear(&point1D);
    linear.calc();
    core_test::expectNear(linear.getPhi(0, 0), 0.75, "linear left value");
    core_test::expectNear(linear.getPhi(0, 1), 0.25, "linear right value");
    core_test::expectNear(linear.getPhi(0, 0) + linear.getPhi(0, 1), 1.0,
                          "linear partition of unity");
    core_test::expectNear(linear.getPhi(1, 0), -0.5, "linear left derivative");
    core_test::expectNear(linear.getPhi(1, 1), 0.5, "linear right derivative");

    mknix::Point pointX(1, 3, 0.5, 0.0, 0.0, 1.0, 2.0);
    addNodes(pointX, nodes);
    mknix::ShapeFunctionLinearX linearX(&pointX);
    linearX.calc();
    core_test::expectNear(linearX.getPhi(0, 0), 0.75, "linear-x left value");
    core_test::expectNear(linearX.getPhi(0, 1), 0.25, "linear-x right value");
}

void testTriangleShapeFunctions()
{
    mknix::Node first(0, 0.0, 0.0, 0.0);
    mknix::Node second(1, 1.0, 0.0, 0.0);
    mknix::Node third(2, 0.0, 1.0, 0.0);
    std::vector<mknix::Node*> nodes = {&first, &second, &third};

    mknix::Point point(2, 4, 0.25, 0.25, 0.0, 0.5, 1.0, 1.0);
    addNodes(point, nodes);
    mknix::ShapeFunctionTriangle triangle(&point);
    triangle.calc();
    core_test::expectNear(triangle.getPhi(0, 0), 0.5, "triangle first barycentric value");
    core_test::expectNear(triangle.getPhi(0, 1), 0.25, "triangle second barycentric value");
    core_test::expectNear(triangle.getPhi(0, 2), 0.25, "triangle third barycentric value");
    core_test::expectNear(triangle.getPhi(0, 0) + triangle.getPhi(0, 1) + triangle.getPhi(0, 2),
                          1.0, "triangle partition of unity");

    mknix::Point point3D(3, 5, 0.25, 0.25, 0.0, 1.0, 1.0, 1.0);
    addNodes(point3D, nodes);
    mknix::ShapeFunctionTriangleSigned signedTriangle(&point3D);
    signedTriangle.calc();
    core_test::expectNear(signedTriangle.getPhi(0, 0), 0.5,
                          "signed triangle first barycentric value");
    core_test::expectNear(signedTriangle.getPhi(0, 1), 0.25,
                          "signed triangle second barycentric value");
    core_test::expectNear(signedTriangle.getPhi(0, 2), 0.25,
                          "signed triangle third barycentric value");
    const double signedSum = signedTriangle.getPhi(0, 0) +
                             signedTriangle.getPhi(0, 1) +
                             signedTriangle.getPhi(0, 2);
    core_test::expectNear(signedSum, 1.0, "signed 3D triangle partition of unity");
}

void testTetrahedronShapeFunction()
{
    mknix::Node first(0, 0.0, 0.0, 0.0);
    mknix::Node second(1, 1.0, 0.0, 0.0);
    mknix::Node third(2, 0.0, 1.0, 0.0);
    mknix::Node fourth(3, 0.0, 0.0, 1.0);
    std::vector<mknix::Node*> nodes = {&first, &second, &third, &fourth};

    mknix::Point point(3, 6, 0.25, 0.25, 0.25, 1.0 / 6.0, 1.0, 1.0);
    addNodes(point, nodes);
    mknix::ShapeFunctionTetrahedron tetrahedron(&point);
    tetrahedron.calc();

    double sum = 0.0;
    for (size_t i = 0; i < 4; ++i)
    {
        core_test::expectNear(tetrahedron.getPhi(0, i), 0.25, "tetrahedron barycentric value");
        sum += tetrahedron.getPhi(0, i);
    }
    core_test::expectNear(sum, 1.0, "tetrahedron partition of unity");
}

void testMeshfreeShapeFunctions()
{
    mknix::Node first(0, 0.0, 0.0, 0.0);
    mknix::Node second(1, 1.0, 0.0, 0.0);
    mknix::Node third(2, 0.0, 1.0, 0.0);
    mknix::Node fourth(3, 0.0, 0.0, 1.0);
    std::vector<mknix::Node*> planarNodes = {&first, &second, &third};
    std::vector<mknix::Node*> spatialNodes = {&first, &second, &third, &fourth};

    mknix::Point mlsPoint(2, 7, 0.2, 0.2, 0.0, 1.0, 2.0);
    addNodes(mlsPoint, planarNodes);
    double alpha2D = 2.0;
    double dc2D = 1.0;
    mknix::ShapeFunctionMLS mls(3, 1, 1, alpha2D, dc2D, &mlsPoint);
    mls.calc();
    core_test::expectNear(mls.getPhi(0, 0) + mls.getPhi(0, 1) + mls.getPhi(0, 2),
                          1.0, "MLS partition of unity");

    mknix::Point mls2DPoint(2, 8, 0.2, 0.2, 0.0, 1.0, 2.0);
    addNodes(mls2DPoint, planarNodes);
    mknix::ShapeFunctionMLS2D mls2D(3, 1, 1, alpha2D, dc2D, &mls2DPoint);
    mls2D.calc();
    core_test::expectNear(mls2D.getPhi(0, 0) + mls2D.getPhi(0, 1) + mls2D.getPhi(0, 2),
                          1.0, "2D MLS partition of unity");

    mknix::Point rbfPoint(2, 9, 0.2, 0.2, 0.0, 1.0, 2.0);
    addNodes(rbfPoint, planarNodes);
    double q = 0.5;
    mknix::ShapeFunctionRBF rbf(3, 0, 0, alpha2D, dc2D, q, &rbfPoint);
    rbf.calc();
    for (size_t i = 0; i < 3; ++i)
    {
        core_test::expect(std::isfinite(rbf.getPhi(0, i)), "RBF shape value is finite");
    }

    mknix::Point mls3DPoint(3, 10, 0.1, 0.1, 0.1, 1.0, 2.0);
    addNodes(mls3DPoint, spatialNodes);
    double alpha3D = 2.0;
    double dc3D = 1.0;
    mknix::ShapeFunctionMLS3D mls3D(4, 1, 1, alpha3D, dc3D, &mls3DPoint);
    mls3D.calc();
    double sum3D = 0.0;
    for (size_t i = 0; i < 4; ++i)
    {
        const double value = mls3D.getPhi(0, i);
        core_test::expect(std::isfinite(value), "3D MLS shape value is finite");
        sum3D += value;
    }
    core_test::expectNear(sum3D, 1.0, "3D MLS partition of unity");
}

}

void testShapeFunctions()
{
    testShapeFunctionBase();
    testLinearShapeFunctions();
    testTriangleShapeFunctions();
    testTetrahedronShapeFunction();
    testMeshfreeShapeFunctions();
}