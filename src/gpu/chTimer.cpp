

#include "chTimer.h"
#include <fstream>
#include <ios>


inline void
chTimerGetTime( chTimerTimestamp *p )
{
    clock_gettime(CLOCK_MONOTONIC, p);
}

inline double
chTimerElapsedTime( chTimerTimestamp *pStart, chTimerTimestamp *pStop )
{ //returns in seconds up to the nanosecond resolution
    return (double) (pStop->tv_sec - pStart->tv_sec) +
                    (pStop->tv_nsec - pStart->tv_nsec)/1e9;
}

void cpuTick(cpuClock *ck){
  chTimerGetTime(&ck->start);
}

void cpuTock(cpuClock *ck, std::string function_name){
  chTimerGetTime(&ck->end);
  double nano = chTimerElapsedTime(&ck->start, &ck->end);
  ck->elapsedMicroseconds = nano * 1e6 ;
  std::cout << function_name <<" Time Consumed: " << ck->elapsedMicroseconds  << " microseconds" <<std::endl;
}

//#endif

inline double
chTimerBandwidth( chTimerTimestamp *pStart, chTimerTimestamp *pStop, double cBytes )
{
    double et = chTimerElapsedTime( pStart, pStop );
    return cBytes / et;
}

void saveClock(cpuClock *ck, std::string file_path ,std::string message){
    std::ofstream s_out(file_path.c_str());
    s_out << message << std::endl;
    s_out << message <<" Time Consumed: " << ck->elapsedMicroseconds  << " microseconds" <<std::endl;
    s_out.close();
}
