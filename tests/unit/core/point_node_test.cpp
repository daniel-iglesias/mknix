#include "core_test_support.h"

#include <core/node.h>
#include <core/point.h>

void testPointAndNode()
{
    mknix::Point point(2, 10, 0.25, 0.5, 0.0, 1.5, 0.5);
    mknix::Point other(11, 1.25, 0.5, 0.0);
    core_test::expect(point.getDim() == 2, "Point stores its spatial dimension");
    core_test::expect(point.getNumber() == 10, "Point stores its identifier");
    core_test::expectNear(point.distance(other), 1.0, "Point distance");

    point.setShapeFunType("MLS");
    point.setDc(0.75);
    point.setDc(0.25);
    core_test::expect(point.getShapeFunType() == "MLS", "Point stores shape-function type");

    mknix::Node first(1, 0.0, 0.0, 0.0);
    mknix::Node second(2, 0.5, 0.0, 0.0);
    mknix::Node outside(3, 1.5, 0.0, 0.0);
    std::vector<mknix::Node*> domain = {&first, &second, &outside};
    mknix::Point supportPoint(2, 20, 0.0, 0.0, 0.0, 1.0, 1.0);
    supportPoint.findSupportNodes(domain);
    core_test::expect(supportPoint.getSupportSize() == 2, "Point selects nodes within its radius");
    core_test::expect(supportPoint.getSupportNodeNumber(0, 0) == 1, "Point exposes support node id");

    first.addWeight(0.25);
    first.addWeight(0.75);
    first.setUx(0.2);
    first.setUy(-0.1);
    first.setUz(0.4);
    first.setqt(315.0);
    first.setThermalNumber(7);
    core_test::expectNear(first.getWeight(), 1.0, "Node accumulates weight");
    core_test::expectNear(first.getConf(0), 0.2, "Node x configuration");
    core_test::expectNear(first.getConf(1), -0.1, "Node y configuration");
    core_test::expectNear(first.getConf(2), 0.4, "Node z configuration");
    core_test::expectNear(first.getTemp(), 315.0, "Node temperature");
    core_test::expectNear(first.getUx(), 0.2, "Node x displacement");
    core_test::expect(first.getThermalNumber() == 7, "Node thermal id");

    mknix::Node copy(first);
    core_test::expect(copy.getNumber() == first.getNumber(), "Node copy preserves id");
    core_test::expectNear(copy.getTemp(), first.getTemp(), "Node copy preserves temperature");
    core_test::expectNear(copy.getConf(0), first.getConf(0), "Node copy preserves configuration");
}