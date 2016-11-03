#ifndef DEVICE_HELPER_H
#define DEVICE_HELPER_H
#include <iostream>
#include <iomanip>
#include <cuda.h>
#include <cuda_runtime.h>
//#include <cuda_runtime_api.h>
//#include <device_launch_parameters.h>
// Define this to turn on error checking
#define CUDA_ERROR_CHECK

//namespace CudaHelper{

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line ){
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",file, line, cudaGetErrorString( err ) );exit( -1 );}
#endif
  return;}
///////////////////////////////////////////////////////////////////////////////
inline void __cudaCheckError( const char *file, const int line ){
#ifdef CUDA_ERROR_CHECK
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err )
    {  fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",file, line, cudaGetErrorString( err ) );  exit( -1 ); }
    // More careful checking. However, this will affect performance.// Comment away if needed.
    /*err = cudaDeviceSynchronize();
    if( cudaSuccess != err ){fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",file, line, cudaGetErrorString( err ) ); exit( -1 );}*/
#endif
  return;}

  /*__device__ double atomicAdd(double* address, double val)
  {
   unsigned long long int* address_as_ull = (unsigned long long int*)address;
   unsigned long long int old = *address_as_ull, assumed;
   do {
     assumed = old;
     old = atomicCAS(address_as_ull, assumed,
     __double_as_longlong(val +__longlong_as_double(assumed)));
   } while (assumed != old);
   return __longlong_as_double(old);
 }*/




/////////////////////////////////////////////////////////////////////////////
void checkGPUMemory()
{ size_t free_byte0 ;
	size_t total_byte0 ;
	if ( cudaSuccess != cudaMemGetInfo( &free_byte0, &total_byte0 ) ){
				std::cout <<"Error: cudaMemGetInfo fails" << std::endl; exit(1);}
	double free_db0 = ((double)free_byte0)/(1024.0*1024.0) ;
	double total_db0 = ((double)total_byte0)/(1024.0*1024.0) ;
	double used_db0 = total_db0 - free_db0 ;
	double memperc = (100 * used_db0)/total_db0;
	std::cout << std::fixed << std::setprecision(1);
	std::cout << std::endl << "---- GPU memory usage : "<< memperc << "% ---" << std::endl;
	std::cout << " used = " << used_db0 <<" MB" << std::endl;
	std::cout << " free = " << free_db0 <<" MB" << std::endl;
	std::cout << " total = " << total_db0 << " MB" << std::endl;
	std::cout << "----------------------------------" << std::endl;}

struct cudaClock
{
  cudaEvent_t start, stop;
  float gpu_time, last_measure;
};

void cudaTick(cudaClock *ck)
{
  cudaEventCreate(&(ck->start));
  cudaEventCreate(&(ck->stop));
  cudaEventRecord(ck->start,0);
}

double cudaTock(cudaClock *ck)
{//modified to suit mknix
  cudaEventRecord(ck->stop, 0);
  cudaEventSynchronize(ck->stop);
  cudaEventElapsedTime(&(ck->gpu_time), ck->start, ck->stop);;
//  std::cout << "GPU clock measured "<<  ck->gpu_time *1000.0f << " microseconds" << std::endl;
  cudaEventDestroy(ck->start); //cleaning up a bit
  cudaEventDestroy(ck->stop);
  ck->last_measure = ck->gpu_time;//not really suing it yet
  ck->gpu_time = 0.0f;
  return ck->last_measure * 1000.0;
}


//}

#endif //DEVICE_HELPER_H
