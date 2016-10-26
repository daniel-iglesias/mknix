#include "chTimer.h"
#include <stdio.h>
#include "gpu_assembly.cuh"


int main(){
	cpuClock ck;
	cpuTick(&ck);
	std::cout << "quick Hello Timer"<< std::endl;
	cpuTock(&ck);
	w_function();
 	return 42;
}
