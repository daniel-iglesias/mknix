#include <simulation/analysisdynamic.h>
#include <simulation/analysisstatic.h>
#include <simulation/analysisthermaldynamic.h>
#include <simulation/analysisthermalstatic.h>
#include <simulation/analysisthermomechanicaldynamic.h>

#include <stdexcept>
#include <string>

void testAnalysisTypes()
{
    mknix::AnalysisStatic staticAnalysis;
    mknix::AnalysisDynamic dynamicAnalysis;
    mknix::AnalysisThermalStatic thermalStaticAnalysis;
    mknix::AnalysisThermalDynamic thermalDynamicAnalysis;
    mknix::AnalysisThermoMechanicalDynamic thermoMechanicalAnalysis;

    if (staticAnalysis.type() != "STATIC" ||
        dynamicAnalysis.type() != "DYNAMIC" ||
        thermalStaticAnalysis.type() != "THERMALSTATIC" ||
        thermalDynamicAnalysis.type() != "THERMAL" ||
        thermoMechanicalAnalysis.type() != "THERMOMECHANICALDYNAMIC")
    {
        throw std::runtime_error("an analysis returned an unexpected type identifier");
    }
}