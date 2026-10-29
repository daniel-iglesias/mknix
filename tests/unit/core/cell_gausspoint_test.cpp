#include "core_test_support.h"

#include <core/cell.h>
#include <core/cellboundary.h>
#include <core/cellboundarylinear.h>
#include <core/cellrect.h>
#include <core/celltetrahedron.h>
#include <core/celltriang.h>
#include <core/elemtetrahedron.h>
#include <core/elemtriangle.h>
#include <core/gausspoint2D.h>
#include <core/gausspoint3D.h>
#include <core/gausspointboundary.h>
#include <core/material.h>
#include <core/node.h>
#include <core/point.h>
#include <core/gausspoint.h>

#ifdef HAVE_VTK
#include <core/compbar.h>
#endif

namespace
{

class InspectableCell : public mknix::Cell
{
public:
    void addPoint(mknix::Point* point)
    {
        bodyPoints.push_back(point);
    }

    void gnuplotOut(std::ofstream&, std::ofstream&) override
    {
    }
};

class InspectableCellRect : public mknix::CellRect
{
public:
    using CellRect::CellRect;

    size_t integrationPointCount() const
    {
        return gPoints.size();
    }

    mknix::GaussPoint* integrationPoint(size_t index)
    {
        return gPoints.at(index);
    }
};

void testCellBase()
{
    mknix::Point first(0, 0.0, 0.0, 0.0);
    mknix::Point second(1, 0.5, 0.0, 0.0);
    mknix::Point third(2, 0.0, 0.5, 0.0);
    InspectableCell cell;
    cell.addPoint(&first);
    cell.addPoint(&second);
    cell.addPoint(&third);
    core_test::expect(cell.getNumberofNodes() == 3, "Cell counts its geometry points");
    const std::vector<int> nodeNumbers = cell.getNodeNumbers();
    core_test::expect(nodeNumbers.size() == 3, "Cell exposes connectivity length");
    core_test::expect(nodeNumbers.at(1) == 1, "Cell exposes connectivity ids");
    mknix::Material replacement;
    core_test::expect(!cell.setMaterialIfLayer(replacement, 0.1),
                      "Cell does not select a layer when point spacing exceeds thickness");
}

void testTriangleCells()
{
    mknix::Material material;
    mknix::Node first(0, 0.0, 0.0, 0.0);
    mknix::Node second(1, 1.0, 0.0, 0.0);
    mknix::Node third(2, 0.0, 1.0, 0.0);

    mknix::CellTriang integrationCell(material, "EFG", 2.0, 1, &first, &second, &third);
    core_test::expect(integrationCell.getNumberofNodes() == 3, "CellTriang has three vertices");
    const std::vector<int> integrationConnectivity = integrationCell.getNodeNumbers();
    core_test::expect(integrationConnectivity.size() == 3, "CellTriang connectivity size");

    mknix::ElemTriangle element(material, 2.0, 1, &first, &second, &third);
    std::vector<mknix::Node*> nodes = {&first, &second, &third};
    element.initialize(nodes);
    element.computeShapeFunctions();
    core_test::expect(element.getNumberofNodes() == 3, "ElemTriangle has three vertices");
    const std::vector<int> elementConnectivity = element.getNodeNumbers();
    core_test::expect(elementConnectivity.at(0) == 0 && elementConnectivity.at(2) == 2,
                      "ElemTriangle preserves node connectivity");
}

void testTetrahedronCells()
{
    mknix::Material material;
    mknix::Node first(0, 0.0, 0.0, 0.0);
    mknix::Node second(1, 1.0, 0.0, 0.0);
    mknix::Node third(2, 0.0, 1.0, 0.0);
    mknix::Node fourth(3, 0.0, 0.0, 1.0);

    mknix::CellTetrahedron integrationCell(material, "EFG", 2.0, 1,
                                           &first, &second, &third, &fourth);
    core_test::expect(integrationCell.getNumberofNodes() == 4,
                      "CellTetrahedron has four vertices");

    mknix::ElemTetrahedron element(material, 2.0, 1,
                                   &first, &second, &third, &fourth);
    std::vector<mknix::Node*> nodes = {&first, &second, &third, &fourth};
    element.initialize(nodes);
    element.computeShapeFunctions();
    core_test::expect(element.getNumberofNodes() == 4, "ElemTetrahedron has four vertices");
    core_test::expect(element.getNodeNumbers().at(3) == 3,
                      "ElemTetrahedron preserves node connectivity");
}

void testRectangularAndBoundaryCells()
{
    mknix::Material material;
    InspectableCellRect rectangle(material, "EFG", 2.0, 1,
                                  0.0, 0.0, 1.0, 0.0,
                                  1.0, 1.0, 0.0, 1.0,
                                  0.5, 0.5, 0.0, 0.0, 1.0, 1.0);
    core_test::expect(rectangle.integrationPointCount() == 1,
                      "CellRect creates one integration point for one-point quadrature");
    core_test::expectNear(rectangle.integrationPoint(0)->getX(), 0.5,
                          "CellRect integration-point x coordinate");
    core_test::expectNear(rectangle.integrationPoint(0)->getY(), 0.5,
                          "CellRect integration-point y coordinate");

    InspectableCellRect hexahedron(material, "EFG", 2.0, 1,
                                   0.0, 0.0, 0.0,
                                   1.0, 0.0, 0.0,
                                   1.0, 1.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 1.0,
                                   1.0, 1.0, 1.0,
                                   0.0, 1.0, 1.0,
                                   0.5, 0.5, 0.5,
                                   0.0, 0.0, 0.0,
                                   1.0, 1.0, 1.0);
    core_test::expect(hexahedron.integrationPointCount() == 1,
                      "3D CellRect creates one integration point");
    core_test::expectNear(hexahedron.integrationPoint(0)->getZ(), 0.5,
                          "3D CellRect integration-point z coordinate");

    mknix::CellBoundary boundary("EFG", 2.0, 1);
    std::vector<mknix::Node*> noNodes;
    boundary.initialize(noNodes);
    boundary.computeShapeFunctions();

    mknix::Node start(0, 0.0, 0.0, 0.0);
    mknix::Node end(1, 1.0, 0.0, 0.0);
    mknix::CellBoundaryLinear boundaryLine("FEM", 2.0, 1, &start, &end);
    std::vector<mknix::Node*> nodes = {&start, &end};
    boundaryLine.initialize(nodes);
    boundaryLine.computeShapeFunctions();
}

void testGaussPointClasses()
{
    mknix::Material material;
    mknix::GaussPoint2D point2D(1.0, 1.0, 0.5, &material, 0, 0.25, 0.25, 1.0, false);
    core_test::expect(point2D.getDim() == 2, "GaussPoint2D dimension");
    core_test::expectNear(point2D.getX(), 0.25, "GaussPoint2D x coordinate");

    mknix::GaussPoint3D point3D(1.0, 1.0, 1.0 / 6.0, &material, 0,
                                0.25, 0.25, 0.25, 1.0, false);
    core_test::expect(point3D.getDim() == 3, "GaussPoint3D dimension");
    core_test::expectNear(point3D.getZ(), 0.25, "GaussPoint3D z coordinate");

    mknix::GaussPointBoundary boundaryPoint(2, 1.0, 1.0, 0.5, 0, 0.25, 1.0);
    core_test::expect(boundaryPoint.getDim() == 2, "GaussPointBoundary dimension");
    core_test::expectNear(boundaryPoint.getX(), 0.25, "GaussPointBoundary x coordinate");

    mknix::Node left(0, 0.0, 0.0, 0.0);
    mknix::Node right(1, 1.0, 0.0, 0.0);
    mknix::Node top(2, 0.0, 1.0, 0.0);
    mknix::ElemTriangle element(material, 1.0, 1, &left, &right, &top);
    std::vector<mknix::Node*> triangleNodes = {&left, &right, &top};
    element.initialize(triangleNodes);
    element.computeShapeFunctions();
}

#ifdef HAVE_VTK
void testCompBar()
{
    mknix::Node first(0, 0.0, 0.0, 0.0);
    mknix::Node second(1, 1.0, 0.0, 0.0);
    mknix::CompBar bar(1, &first, &second);
    bar.updatePoints();
}
#endif

}

void testCellsAndGaussPoints()
{
    testCellBase();
    testTriangleCells();
    testTetrahedronCells();
    testRectangularAndBoundaryCells();
    testGaussPointClasses();
#ifdef HAVE_VTK
    testCompBar();
#endif
}