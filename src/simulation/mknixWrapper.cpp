#include "mknixWrapper.h"

#include "simulation.h"

struct mknix::MknixWrapper::SimulationWrapper {
  mknix::Simulation theSimulation;
};

mknix::MknixWrapper::MknixWrapper(const std::string& input_file_name, const int input_detail)
    : _sim(new mknix::MknixWrapper::SimulationWrapper()) {
  config(input_file_name, input_detail);
}

mknix::MknixWrapper::~MknixWrapper() = default;

void mknix::MknixWrapper::config(const std::string& input_file_name, int output_detail) {
  _sim->theSimulation.setOutputFilesDetail(output_detail);
  _sim->theSimulation.inputFromFile(input_file_name);
};

void mknix::MknixWrapper::init(double temperatures) {
  _sim->theSimulation.setInitialTemperatures(temperatures);
  _sim->theSimulation.init();
}

void mknix::MknixWrapper::run(double* heatFluence, double* temperatures) {
  _sim->theSimulation.solveStep(heatFluence, temperatures);
}