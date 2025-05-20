#pragma once

#include <memory>

namespace mknix {

class MknixWrapper {
 public:
  MknixWrapper(const std::string& input_file_name, const int input_detail = 0);
  virtual ~MknixWrapper();
  
  void init(double temperatures);
  
  void run(double* heatFluence, double* temperatures);
 
 private:
  void config(const std::string& input_file_name, int output_detail);
  
  struct SimulationWrapper;
  std::unique_ptr<SimulationWrapper> _sim;
};

}  // namespace mknix
